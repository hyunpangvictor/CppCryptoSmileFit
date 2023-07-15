#ifndef QF633_CODE_MSG_H
#define QF633_CODE_MSG_H

#include <cstdint>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include "Date.h"
#include <sstream>

enum class Moneyness
{
    ITM,
    ATM,
    OTM
};

struct TickData
{
    std::string ContractName;
    double BestBidPrice;
    double BestBidAmount;
    double BestBidIV;
    double BestAskPrice;
    double BestAskAmount;
    double BestAskIV;
    double MarkPrice;
    double MarkIV;
    std::string UnderlyingIndex;
    double UnderlyingPrice;
    double LastPrice;
    double OpenInterest;
    uint64_t LastUpdateTimeStamp;
    std::string Underlying;
    datetime_t ExpiryDate;
    double Strike;
    std::string OptionType;
    Moneyness moneyness;

    // defining constructor for to reduce the copy and paste of objects when pushing into array
    TickData(
        std::string ContractName,
        double BestBidPrice,
        double BestBidAmount,
        double BestBidIV,
        double BestAskPrice,
        double BestAskAmount,
        double BestAskIV,
        double MarkPrice,
        double MarkIV,
        std::string UnderlyingIndex,
        double UnderlyingPrice,
        double LastPrice,
        double OpenInterest,
        uint64_t LastUpdateTimeStamp)
        : ContractName(ContractName),
          BestBidPrice(BestBidPrice),
          BestBidAmount(BestBidAmount),
          BestBidIV(BestBidIV/100),
          BestAskPrice(BestAskPrice),
          BestAskAmount(BestAskAmount),
          BestAskIV(BestAskIV/100),
          MarkPrice(MarkPrice),
          MarkIV(MarkIV/100),
          UnderlyingIndex(UnderlyingIndex),
          UnderlyingPrice(UnderlyingPrice),
          LastPrice(LastPrice),
          OpenInterest(OpenInterest),
          LastUpdateTimeStamp(LastUpdateTimeStamp),
          Underlying(""),
          ExpiryDate(),
          Strike(0.0),
          OptionType(""),
          moneyness(Moneyness::ATM) // Initialize moneyness to ATM
    {
        ExtractContractInfo();
        // logging it out just to ensure the data is being added correctly
        // LogTickData();
    }

    void LogMoneyness() const
    {
        std::string moneynessStr;

        switch (moneyness)
        {
        case Moneyness::ITM:
            moneynessStr = "ITM";
            break;
        case Moneyness::ATM:
            moneynessStr = "ATM";
            break;
        case Moneyness::OTM:
            moneynessStr = "OTM";
            break;
        default:
            moneynessStr = "ATM";
            break;
        }

        std::cout << "Moneyness: " << moneynessStr << std::endl;
    }
    void LogTickData() const
    {
        // logging it out just to ensure the data is being added correctly
        std::cout << "-------------" << std::endl;
        std::cout << "Tick Data:" << std::endl;
        std::cout << "-------------" << std::endl;
        std::cout << "ContractName : " << ContractName << std::endl;
        std::cout << "BestBidPrice : " << BestBidPrice << std::endl;
        std::cout << "BestBidAmount : " << BestBidAmount << std::endl;
        std::cout << "BestBidIV : " << BestBidIV << std::endl;
        std::cout << "BestAskPrice : " << BestAskPrice << std::endl;
        std::cout << "BestAskAmount : " << BestAskAmount << std::endl;
        std::cout << "BestAskIV : " << BestAskIV << std::endl;
        std::cout << "MarkPrice : " << MarkPrice << std::endl;
        std::cout << "MarkIV : " << MarkIV << std::endl;
        std::cout << "UnderlyingIndex : " << UnderlyingIndex << std::endl;
        std::cout << "UnderlyingPrice : " << UnderlyingPrice << std::endl;
        std::cout << "LastPrice : " << LastPrice << std::endl;
        std::cout << "OpenInterest : " << OpenInterest << std::endl;
        std::cout << "LastUpdateTimeStamp : " << LastUpdateTimeStamp << std::endl;
        std::cout << "Underlying : " << Underlying << std::endl;
        std::cout << "ExpiryDate : " << ExpiryDate << std::endl;
        std::cout << "Strike : " << Strike << std::endl;
        std::cout << "OptionType : " << OptionType << std::endl;
        LogMoneyness();
        std::cout << "-------------" << std::endl;
    }

    void ExtractContractInfo()
    {
        size_t firstDash = ContractName.find('-');
        size_t secondDash = ContractName.find('-', firstDash + 1);
        size_t thirdDash = ContractName.find('-', secondDash + 1);

        if (firstDash != std::string::npos && secondDash != std::string::npos && thirdDash != std::string::npos)
        {
            Underlying = ContractName.substr(0, firstDash);
            ExpiryDate = GetExpiry(ContractName.substr(firstDash + 1, secondDash - firstDash - 1));
            Strike = std::stod(ContractName.substr(secondDash + 1, thirdDash - secondDash - 1));
            OptionType = ContractName.substr(thirdDash + 1);
        }

        if (OptionType == "C")
        {
            if (Strike < UnderlyingPrice)
            {
                moneyness = Moneyness::ITM;
            }
            else if (Strike == UnderlyingPrice)
            {
                moneyness = Moneyness::ATM;
            }
            else
            {
                moneyness = Moneyness::OTM;
            }
        }
        else if (OptionType == "P")
        {
            if (Strike > UnderlyingPrice)
            {
                moneyness = Moneyness::ITM;
            }
            else if (Strike == UnderlyingPrice)
            {
                moneyness = Moneyness::ATM;
            }
            else
            {
                moneyness = Moneyness::OTM;
            }
        }
    }

    datetime_t GetExpiry(const std::string &str_datetime)
    {

        std::istringstream iss(str_datetime);
        std::tm time = {};
        iss >> std::get_time(&time, "%d%b%y");

        datetime_t datetime;
        datetime.year = time.tm_year + 1900; // tm_year is years since 1900
        datetime.month = time.tm_mon + 1;    // tm_mon is 0-based, so add 1
        datetime.day = time.tm_mday;

        return datetime;
    }
};

struct Msg
{
    uint64_t timestamp{};
    bool isSnap;
    bool isSet = false;
    std::vector<TickData> Updates;

    // Log function to print Msg object
    void LogMsg() const
    {
        // Iterate over the members of the Msg struct
        std::cout << "-------------" << std::endl;
        std::cout << "Msg object:" << std::endl;
        std::cout << "-------------" << std::endl;
        std::cout << "timestamp: " << timestamp << std::endl;
        std::cout << "isSnap: " << isSnap << std::endl;
        std::cout << "isSet: " << isSet << std::endl;
        std::cout << "Updates size: " << Updates.size() << std::endl;
        std::cout << "-------------" << std::endl;
    }
};

#endif // QF633_CODE_MSG_H
