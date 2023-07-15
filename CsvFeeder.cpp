#include "CsvFeeder.h"
#include "date/date.h"
#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>

uint64_t TimeToUnixMS(std::string ts)
{
  std::istringstream in{ts};
  std::chrono::system_clock::time_point tp;
  in >> date::parse("%FT%T", tp);
  const auto timestamp =
      std::chrono::time_point_cast<std::chrono::milliseconds>(tp)
          .time_since_epoch()
          .count();
  return timestamp;
}

double ConvertToDouble(const std::string &str)
{
  try
  {
    return std::stod(str);
  }
  catch (...)
  {
    // if fields are empty set to nan
    return std::numeric_limits<double>::quiet_NaN();
  }
}

bool ReadHeader(std::ifstream &file, Msg &msg, std::map<std::string, int> &column_pos)
{
  // reading the file for the first time
  std::string line;
  std::stringstream ss(line);
  std::string token;

  // read the first line and set column positions
  std::getline(file, line);
  ss.str(line);

  // std::cout << "reading header position:" << file.tellg() << std::endl;
  int pos = 0;
  while (std::getline(ss, token, ','))
  {
    // std::cout << "assigning header token : " << token << " in position : " << pos << std::endl;
    column_pos[token] = pos;
    pos++;
  }

  // read the next line and insert into row vector
  std::getline(file, line);

  std::stringstream ssheader(line);
  std::vector<std::string> row_vec;
  std::string token2;
  // std::cout << "current position:" << file.tellg() << std::endl;
  // pull out row data from csv
  while (std::getline(ssheader, token2, ','))
  {
    row_vec.emplace_back(std::move(token2));
  }

  // initialise the first time stamp
  msg.timestamp = TimeToUnixMS(row_vec[column_pos["time"]]);

  file.seekg(-(line.size() + 1), std::ios::cur);
  std::cout << "read finished first time stamp reset current position to ->" << file.tellg() << std::endl;

  return true;
}

void ReplaceTickData(Msg &msg, const std::string &contractName, const TickData &newTickData)
{
  for (auto &tickData : msg.Updates)
  {
    if (tickData.ContractName == contractName)
    {
      tickData = newTickData;
      break; // Found the contract, no need to continue searching
    }
  }
}

bool ReadNextMsg(std::ifstream &file, Msg &msg, std::map<std::string, int> &column_pos)
{
  // TODO: your implementation to read file and create the next Msg into the
  // std::cout << "[ReadNextMsg] start..." << std::endl;
  // Get the initial position of the file pointer
  std::streampos initialPosition = file.tellg();

  if (file.eof())
  {
    return false;
  }
  else if (msg.isSet == false)
  {
    // read the file for the first time
    msg.isSet = true;
    // std::cout << "reading header columns..." << std::endl;
    return ReadHeader(file, msg, column_pos);
  }

  std::string line;
  uint64_t current_ts = 0U; // currently processing row's ts
  uint64_t last_ts = 0U;    // previously processed row's ts
  uint64_t first_snap_ts = 0U;

  // Print the positions
  // std::cout << "reading row data Position before while loop: " << initialPosition << std::endl;
  while (current_ts == last_ts || msg.Updates.empty())
  {

    // read the next line and insert into row vector
    std::getline(file, line);
    if (line.empty())
    {
      return false;
    }

    std::stringstream ss2(line);
    std::vector<std::string> row_vec;
    std::string token;

    // std::cout << "current position:" << file.tellg() << std::endl;

    // pull out row data from csv
    while (std::getline(ss2, token, ','))
    {
      row_vec.emplace_back(std::move(token));
    }

    // check if timestamp is the same as last read
    msg.timestamp = TimeToUnixMS(row_vec[column_pos["time"]]);
    current_ts = msg.timestamp;
    // std::cout << "updated current ts with row data: " << current_ts << std::endl;

    if (msg.Updates.empty())
    {
      if (row_vec[column_pos["msgType"]] != "snap")
      {
        // skip rows until first snap message found
        continue;
      }
      else
      {
        first_snap_ts = current_ts;
      }
    }

    if (current_ts != last_ts && last_ts != 0U)
    {
      // time changed we need to step out
      // so we can re process this new time stamp
      file.seekg(-(line.size() + 1), std::ios::cur);
      // std::cout << "updated last ts: " << last_ts << std::endl;
      // std::cout << "change in time stamp reset current position to ->" << file.tellg() << std::endl;
      break;
    }

    // update msg
    // add tickdata to msg
    msg.isSnap = row_vec[column_pos["msgType"]] == "snap";
    if (msg.isSnap)
    {
      // std::cout << "msg is a snap proceed to check if its a subsequent snaps which require clear..." << std::endl;
      if (first_snap_ts != TimeToUnixMS(row_vec[column_pos["time"]]))
      {
        // std::cout << "first snap time is : " << first_snap_ts << " different from current row time proceed to clear... " << std::endl;
        msg.Updates.clear();
      }
      // add tickdata to msg after droping all previous snaps
      // std::cout << "adding tick data updates..." << std::endl;
      msg.Updates.emplace_back(
          row_vec[column_pos["contractName"]],
          ConvertToDouble(row_vec[column_pos["bestBid"]]),
          ConvertToDouble(row_vec[column_pos["bestBidAmount"]]),
          ConvertToDouble(row_vec[column_pos["bestBidIV"]]),
          ConvertToDouble(row_vec[column_pos["bestAsk"]]),
          ConvertToDouble(row_vec[column_pos["bestAskAmount"]]),
          ConvertToDouble(row_vec[column_pos["bestAskIV"]]),
          ConvertToDouble(row_vec[column_pos["markPrice"]]),
          ConvertToDouble(row_vec[column_pos["markIV"]]),
          row_vec[column_pos["underlyingIndex"]],
          ConvertToDouble(row_vec[column_pos["underlyingPrice"]]),
          ConvertToDouble(row_vec[column_pos["lastPrice"]]),
          ConvertToDouble(row_vec[column_pos["open_interest"]]),
          TimeToUnixMS(row_vec[column_pos["time"]]));
    }
    else
    {
      // std::cout << "msg is not a snap, proceed with replace tick data..." << std::endl;
      ReplaceTickData(msg,
                      row_vec[column_pos["contractName"]],
                      TickData(row_vec[column_pos["contractName"]],
                               ConvertToDouble(row_vec[column_pos["bestBid"]]),
                               ConvertToDouble(row_vec[column_pos["bestBidAmount"]]),
                               ConvertToDouble(row_vec[column_pos["bestBidIV"]]),
                               ConvertToDouble(row_vec[column_pos["bestAsk"]]),
                               ConvertToDouble(row_vec[column_pos["bestAskAmount"]]),
                               ConvertToDouble(row_vec[column_pos["bestAskIV"]]),
                               ConvertToDouble(row_vec[column_pos["markPrice"]]),
                               ConvertToDouble(row_vec[column_pos["markIV"]]),
                               row_vec[column_pos["underlyingIndex"]],
                               ConvertToDouble(row_vec[column_pos["underlyingPrice"]]),
                               ConvertToDouble(row_vec[column_pos["lastPrice"]]),
                               ConvertToDouble(row_vec[column_pos["open_interest"]]),
                               TimeToUnixMS(row_vec[column_pos["time"]])));
    }

    last_ts = msg.timestamp;

    // std::cout << "after adding tick data to updates" << std::endl;
    // std::cout << "last ts: " << last_ts << std::endl;
    // std::cout << "current ts: " << current_ts << std::endl;
    // msg.LogMsg();
  }
  return true;
}

CsvFeeder::CsvFeeder(const std::string ticker_filename,
                     FeedListener feed_listener, std::chrono::minutes interval,
                     TimerListener timer_listener)
    : ticker_file_(ticker_filename), feed_listener_(feed_listener),
      interval_(interval), timer_listener_(timer_listener)
{
  // initialize member variables with input information, prepare for Step()
  // processing

  // std::cout << "File path: "
  //           << std::filesystem::absolute(ticker_filename).string() << std::endl;

  // std::cout << "Is file open: " << ticker_file_.is_open() << std::endl;

  ReadNextMsg(ticker_file_, msg_, column_pos_);
  if (msg_.isSet)
  {
    // initialize interval timer now_ms_
    now_ms_ = msg_.timestamp;
  }
  else
  {
    throw std::invalid_argument("empty message at initialization");
  }
}

bool CsvFeeder::Step()
{
  // std::cout << "[CsvFeeder::Step()] start ..." << std::endl;
  if (msg_.isSet)
  {
    // call feed_listener with the loaded Msg
    feed_listener_(msg_);

    // if current message's timestamp is crossing the given interval, call
    // time_listener, change now_ms_ to the next interval cutoff
    if (now_ms_ < msg_.timestamp)
    {
      timer_listener_(now_ms_);
      now_ms_ += interval_.count();
    }
    // load tick data into Msg
    // if there is no more message from the csv file, return false, otherwise
    // true
    return ReadNextMsg(ticker_file_, msg_, column_pos_);
  }

  return false;
}

CsvFeeder::~CsvFeeder()
{
  // release resource allocated in constructor, if any
}