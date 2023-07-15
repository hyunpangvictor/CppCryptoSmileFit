#include <iostream>
#include <ctime>

#include "CsvFeeder.h"
#include "Msg.h"
#include "VolSurfBuilder.h"
#include "CubicSmile.h"
#include "helper.h"

// create header in output file
void init_output_file(std::ofstream &file)
{
    std::vector<std::string> col_vec{"TIME", "EXPIRY", "FUT_PRICE", "ATM", "BF25", "RR25", "BF10", "RR10", "ERROR",
                                     "BEF_ATM", "BEF_BF25", "BEF_RR25", "BEF_BF10", "BEF_RR10", "BEF_ERROR"};
    // Send column names to the stream
    for (auto j = 0U; j < col_vec.size(); ++j)
    {
        file << col_vec.at(j);
        if (j != col_vec.size() - 1)
            file << ","; // No comma at end of line
    }
    file << "\n";
}

void output_file(std::ofstream &file, const std::map<datetime_t, std::pair<CubicSmile, double>> &smile, const std::string &current_time)
{
    for (const auto &[datetime, smilePair] : smile)
    {
        const CubicSmile &smileObj = smilePair.first;
        const double fittingError = smilePair.second;

        const std::vector<std::pair<double, double>> &strikeMarks = smileObj.GetStrikeMarks();

        // Convert expiry to formatted string
        char dateStrBuffer[30];
        std::sprintf(dateStrBuffer, "%02d-%s-%04d",
             datetime.day, getMonthAbbreviation(datetime.month).c_str(), datetime.year);

        // convert to decimal points
        double v_qd90 = strikeMarks[0].second;
        double v_qd75 = strikeMarks[1].second;
        double atmvol = strikeMarks[2].second;
        double v_qd25 = strikeMarks[3].second;
        double v_qd10 = strikeMarks[4].second;

        double bf25 = ((v_qd75 + v_qd25) / 2.0) - atmvol;
        double bf10 = ((v_qd90 + v_qd10) / 2.0) - atmvol;
        double rr25 = v_qd25 - v_qd75;
        double rr10 = v_qd10 - v_qd90;

        file << current_time << ","
             << std::string(dateStrBuffer) << ","
             << smileObj.future_price << ","
             << atmvol << ","
             << bf25 << ","
             << rr25 << ","
             << bf10 << ","
             << rr10 << ","
             << fittingError << ","
             << smileObj.primer_guess[0] << ","
             << smileObj.primer_guess[1] << ","
             << smileObj.primer_guess[2] << ","
             << smileObj.primer_guess[3] << ","
             << smileObj.primer_guess[4] << ","
             << smileObj.primer_error << "\n";
    }
}

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "Usage: "
                  << argv[0] << " tick_data.csv"
                  << " outputFile.csv" << std::endl;
        return 1;
    }
    const char *ticker_filename = argv[1];
    const char *output_filename = argv[2];

    std::ofstream outputFile(output_filename);
    init_output_file(outputFile);

    VolSurfBuilder<CubicSmile> volBuilder;
    auto feeder_listener = [&volBuilder](const Msg &msg)
    {
        if (msg.isSet)
        {
            volBuilder.Process(msg);
        }
    };

    auto timer_listener = [&](uint64_t now_ms)
    {
        std::string now_time_str = convert_msec_to_utc_string(now_ms);
        // fit smile
        auto smiles = volBuilder.FitSmiles();
        // TODO: stream the smiles and their fitting error to outputFile.csv
        output_file(outputFile, smiles, now_time_str);
    };

    const auto interval = std::chrono::minutes(1); // we call timer_listener at 1 minute interval
    CsvFeeder csv_feeder(ticker_filename,
                         feeder_listener,
                         interval,
                         timer_listener);

    while (csv_feeder.Step())
    {
    }
    return 0;
}