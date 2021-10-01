#ifndef POLAR_H
#define POLAR_H

//headers
#include<vector>
#include<math.h>
#include<assert.h>
#include<string>
#include<algorithm>
#include<numeric>

// class definition
class POLAR {
private:
    //member variables
    unsigned m_N; // codeword length
    unsigned m_M; // mother codeword length
    unsigned m_K; // information length
    std::vector<unsigned> m_F;		// frozen bit positions (including puncture postions)
    std::vector<unsigned> m_P;		// puncture pattern
    std::vector<unsigned> m_I;		// information bit positions
    std::vector<unsigned> m_Q;		// reliability sequence

public://public member functions

    //constructor
    POLAR(uint32_t info_length,
          uint32_t code_length,
          const std::string construct_method = "Huawei Approx",
          const double design_para = 0);
    //encoder
    std::vector<bool> encoder(std::vector<bool>* msg);
    //SC decoder
    std::vector<bool> sc_decoder(std::vector<double>* llr);
    //SCL decoder
    std::vector<bool> scl_decoder(std::vector<double>* llr, std::vector<bool> crc_g, unsigned nL);
    //rate matching
    std::vector<bool> rate_matching(std::vector<bool>* in);
    //rate recovery
    std::vector<double> rate_recovery(std::vector<double>* in);

    //static member functions can be used outside without an object
    static std::vector<double> channel_polarization_huawei_approx(unsigned N);
    template<typename T>
    static std::vector<size_t> sort_indexes(const std::vector<T>& v);

    //crc related
    static std::vector<bool> crc_gen(std::vector<bool>* msg, std::vector<bool> crc_g);
    static bool crc_check_sum(std::vector<bool>* msg, std::vector<bool> crc_g);
    static void crc_division(std::vector<int>* pad_msg, std::vector<int> gen);
    static std::vector<bool> crc_generator(const std::string crc_type);

private: //private member functions

    //encode
    std::vector<bool> polar_encode(std::vector<bool> u);

    //SC decode
    std::vector<bool> sc_node_operations(std::vector<double>* alpha,
                                         std::vector<bool> F,
                                         std::vector<bool>* u_cap);

    //SCL decode
    void scl_node_operations(std::vector<std::vector<double> >* llr,
                             std::vector<std::vector<bool> >* hard_decision,
                             std::vector<bool> F,
                             std::vector<double>* PM,
                             std::vector<std::vector<bool> >* u_cap);
    // core function of SCL: sort and prune sc decoders
    template <typename T>
    inline void sort_and_prune(std::vector<std::vector<T> >* V, std::vector<uint8_t> P);
};

#endif // !POLAR_H
