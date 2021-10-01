// PolarCodeCPPStdClass.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "polar.h"
#include <random>
#include <cstdarg> //for format

using namespace std;

// function declaration
string format(const char* fmt, ...);
vector<double> linspace(double start, double ed, int num);
void example_run_SCL();
void example_run_SC();

int main()
{
    std::cout << "Hello World!\n";
    example_run_SCL();
    //example_run_SC();
}

void example_run_SCL()
{
    unsigned nL = 4;
    int M = 4000; // code word length
    double rate = 1.0 / 3; // code rate

    unsigned K = ceil(M * rate); // information length

    // crc generator
    vector<bool> crc_g = { 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1 };
    unsigned crc_len = crc_g.size() - 1;

    // information length excluding crc
    unsigned A = K - crc_len;

    // instantiate a POLAR object
    POLAR polar = POLAR(K, M);

    // random engines
    default_random_engine random_engine;
    bernoulli_distribution  bern_dist;
    normal_distribution<double> norm_dist(0, 1);

    // Running parameters
    vector<double> EsN0_dB = linspace(-5, -2, 7);
    vector<double> N0;
    for (auto e : EsN0_dB)
        N0.push_back(pow(10.0, -e / 10));

    vector<double> ber(N0.size(), 0), bler(N0.size(), 0);
    vector<unsigned> n_bit_errs(N0.size(), 0), n_blk_errs(N0.size(), 0);

    unsigned n_max_blks = 10000;

    // loop each SNR
    for (unsigned i = 0; i < N0.size(); i++) {
        //print progress
        string str = format("\nNow running EsN0: %.2f dB [%d of %d]", EsN0_dB[i], i + 1, N0.size());
        cout << str << endl;
        unsigned print_len = 0;

        unsigned n_blks_done = 0;
        clock_t tStart = clock(); // timing

        while ((n_blks_done < n_max_blks) && (n_blk_errs[i] < 100)) {
            // generate random bit stream
            vector<bool> msg; msg.reserve(A);
            for (unsigned j = 0; j < A; j++)
                msg.push_back(bern_dist(random_engine));

            // generate CRC and attach it to the message
            vector<bool> crc = POLAR::crc_gen(&msg, crc_g);
            vector<bool> crc_msg = msg; crc_msg.reserve(K);
            crc_msg.insert(crc_msg.end(), crc.begin(), crc.end());

            // polar encoding
            vector<bool> enc = polar.encoder(&crc_msg);

            //rate matching
            vector<bool> rm_enc = polar.rate_matching(&enc);

            // BPSK+ AWGN
            vector<double> r; r.reserve(M);
            for (auto e : rm_enc)
                r.push_back(1 - 2.0 * e + sqrt(N0[i] / 2.0) * norm_dist(random_engine));

            // compute soft bits as LLR
            vector<double> llr; llr.reserve(M);
            for (auto e : r)
                llr.push_back(4.0 * e / N0[i]);

            // rate recovery
            vector<double> rr_llr = polar.rate_recovery(&llr);

            // scl decoding
            vector<bool> msg_cap = polar.scl_decoder(&rr_llr, crc_g, nL);

            // count errors
            unsigned n_errs = 0;
            for (int j = 0; j < A; j++) {
                if (msg[j] != msg_cap[j])
                    n_errs++;
            }

            if (n_errs) {
                n_bit_errs[i] += n_errs;
                n_blk_errs[i]++;
            }

            n_blks_done += 1;

            ber[i] = n_bit_errs[i] * 1.0 / (A * n_blks_done);
            bler[i] = n_blk_errs[i] * 1.0 / n_blks_done;

            // print progress for every 10 blocks
            if (n_blks_done % 10 == 0) {
                str = format("Elapsed time: %.1f seconds, # tx blocks: %d,# error blocks:%d, ber: %.5f, bler %.5f", (clock() - tStart) / (double)CLOCKS_PER_SEC, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
                cout << std::string(print_len, '\b');
                cout << str << flush;
                print_len = str.length();
            }
        }

        // print  progress when one SNR is finished
        str = format("Elapsed time: %.1f seconds, # tx blocks: %d,# error blocks:%d, ber: %.5f, bler %.5f", (clock() - tStart) / (double)CLOCKS_PER_SEC, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
        cout << std::string(print_len, '\b') << str<< flush;
    }

    // print simulation result
    cout << endl;
    cout << "Modulation:" << "BPSK" << endl;
    cout << "[M,R] = [ " << M << "," << rate << "]" << endl;
    cout << "EsN0_dB = [";
    for (auto e : EsN0_dB)
        cout << e << " ";
    cout << "]" << endl;

    cout << "BER = [";
    for (auto e : ber)
        cout << e << " ";
    cout << "]" << endl;

    cout << "BLER = [";
    for (auto e : bler)
        cout << e << " ";
    cout << "]" << endl;
}
void example_run_SC()
{
    unsigned nL = 1; // number of decoders
    int M = 4000; // code word length
    double rate = 1.0 / 3; // code rate

    unsigned K = ceil(M * rate); // information length

    // crc generator
    vector<bool> crc_g = { 1 };
    unsigned crc_len = crc_g.size() - 1;

    // information length excluding crc
    unsigned A = K - crc_len;

    // instantiate a POLAR object
    POLAR polar = POLAR(K, M);

    // random engines
    default_random_engine random_engine;
    bernoulli_distribution  bern_dist;
    normal_distribution<double> norm_dist(0, 1);

    // Running parameters
    vector<double> EsN0_dB = linspace(-5, -2, 7);
    vector<double> N0;
    for (auto e : EsN0_dB)
        N0.push_back(pow(10.0, -e / 10));

    vector<double> ber(N0.size(), 0), bler(N0.size(), 0);
    vector<unsigned> n_bit_errs(N0.size(), 0), n_blk_errs(N0.size(), 0);

    unsigned n_max_blks = 10000;

    // loop each SNR
    for (unsigned i = 0; i < N0.size(); i++) {
        //print progress
        string str = format("\nNow running EsN0: %.2f dB [%d of %d]", EsN0_dB[i], i + 1, N0.size());
        cout << str << endl;
        unsigned print_len = 0;

        unsigned n_blks_done = 0;
        clock_t tStart = clock(); // timing

        while ((n_blks_done < n_max_blks) && (n_blk_errs[i] < 100)) {
            // generate random bit stream
            vector<bool> msg; msg.reserve(A);
            for (unsigned j = 0; j < A; j++)
                msg.push_back(bern_dist(random_engine));

            // generate CRC and attach it to the message
            vector<bool> crc = POLAR::crc_gen(&msg, crc_g);
            vector<bool> crc_msg = msg; crc_msg.reserve(K);
            crc_msg.insert(crc_msg.end(), crc.begin(), crc.end());

            // polar encoding
            vector<bool> enc = polar.encoder(&crc_msg);

            //rate matching
            vector<bool> rm_enc = polar.rate_matching(&enc);

            // BPSK+ AWGN
            vector<double> r; r.reserve(M);
            for (auto e : rm_enc)
                r.push_back(1 - 2.0 * e + sqrt(N0[i] / 2.0) * norm_dist(random_engine));

            // compute soft bits as LLR
            vector<double> llr; llr.reserve(M);
            for (auto e : r)
                llr.push_back(4.0 * e / N0[i]);

            // rate recovery
            vector<double> rr_llr = polar.rate_recovery(&llr);

            // SC decoder
            vector<bool> msg_cap = polar.sc_decoder(&rr_llr);

            // scl decoding
            //vector<bool> msg_cap = polar.scl_decoder(&rr_llr, crc_g, nL);

            // count errors
            unsigned n_errs = 0;
            for (int j = 0; j < A; j++) {
                if (msg[j] != msg_cap[j])
                    n_errs++;
            }

            if (n_errs) {
                n_bit_errs[i] += n_errs;
                n_blk_errs[i]++;
            }

            n_blks_done += 1;

            ber[i] = n_bit_errs[i] * 1.0 / (A * n_blks_done);
            bler[i] = n_blk_errs[i] * 1.0 / n_blks_done;

            // print progress for every 10 blocks
            if (n_blks_done % 10 == 0) {
                str = format("Elapsed time: %.1f seconds, # tx blocks: %d,# error blocks:%d, ber: %.5f, bler %.5f", (clock() - tStart) / (double)CLOCKS_PER_SEC, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
                cout << std::string(print_len, '\b');
                cout << str << flush;
                print_len = str.length();
            }
        }

        // print  progress when one SNR is finished
        str = format("Elapsed time: %.1f seconds, # tx blocks: %d,# error blocks:%d, ber: %.5f, bler %.5f", (clock() - tStart) / (double)CLOCKS_PER_SEC, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
        cout << std::string(print_len, '\b') << str<< flush;
    }

    // print simulation result
    cout << endl;
    cout << "Modulation:" << "BPSK" << endl;
    cout << "[M,R] = [ " << M << "," << rate << "]" << endl;
    cout << "EsN0_dB = [";
    for (auto e : EsN0_dB)
        cout << e << " ";
    cout << "]" << endl;

    cout << "BER = [";
    for (auto e : ber)
        cout << e << " ";
    cout << "]" << endl;

    cout << "BLER = [";
    for (auto e : bler)
        cout << e << " ";
    cout << "]" << endl;
}
string format(const char* fmt, ...) {
    int size = 512;
    char* buffer = 0;
    buffer = new char[size];
    va_list vl;
    va_start(vl, fmt);
    int nsize = vsnprintf(buffer, size, fmt, vl);
    if (size <= nsize) { //fail delete buffer and try again
        delete[] buffer;
        buffer = 0;
        buffer = new char[nsize + 1]; //+1 for /0
        nsize = vsnprintf(buffer, size, fmt, vl);
    }
    string ret(buffer);
    va_end(vl);
    delete[] buffer;
    return ret;
}
vector<double> linspace(double start, double ed, int num) {
    // catch rarely, throw often
    if (num < 2) {
        throw new exception();
    }
    int partitions = num - 1;
    vector<double> pts;
    // length of each segment
    double length = (ed - start) / partitions;
    // first, not to change
    pts.push_back(start);
    for (int i = 1; i < num - 1; i++) {
        pts.push_back(start + i * length);
    }
    // last, not to change
    pts.push_back(ed);
    return pts;
}
