#include "polar.h"

using std::string;
using std::vector;

// public number functions------------------------------------------------------------------------------------------------
POLAR::POLAR(uint32_t info_length, uint32_t code_length, const string construct_method, const double design_para)
{
    //---------------------------------------------------------------------------------
    // polar construction
    // Input:
    //		info_length: information bits length
    //		code_length: coded bit length
    // 		construct_method: "Huawei approx", "BP", "GA"
    // 		design_para: N/A for"Huawei approx"; erasure probability for "BP"; design SNR in dB for "GA"
    // Dr. J Mao 2021 Sep
    //---------------------------------------------------------------------------------

    m_M = code_length;					// codeword length
    m_N = static_cast<unsigned>(pow(2, ceil(log2(code_length)))); // mother codeword length
    m_K = info_length;					// information length

    // obtain sub-channel's reliability
    vector<double> W(m_N, 0);
    if (construct_method.compare("Huawei Approx") == 0) {
        W = channel_polarization_huawei_approx(m_N);
    }
    else {
        //other methods to be added
        perror("unknown construction methods!");
    }

    // reliability sequence
    for (auto e : sort_indexes(W)) {
        m_Q.push_back(e);
    }

    // puncture pattern
    m_P.reserve(m_N - m_M);
    m_P = vector<unsigned>(m_Q.begin(), m_Q.begin() + (m_N - m_M));
    sort(m_P.begin(), m_P.end());

    //frozen bits positions
    m_F.reserve(m_N - m_K);
    m_F = vector<unsigned>(m_Q.begin(), m_Q.begin() + (m_N - m_K));
    sort(m_F.begin(), m_F.end());

    // information bit positions
    m_I.reserve(m_K);
    m_I = vector<unsigned>(m_Q.end() - m_K, m_Q.end());
    sort(m_I.begin(), m_I.end());
}
vector<bool> POLAR::encoder(vector<bool>* msg)
{
    //---------------------------------------------------------------------------------
    // Polar code encoder
    // Input:
    //		msg: a bvec of size K
    // Return:
    //		the polar code encoded message
    // Dr. J Mao 2021 Sep
    //---------------------------------------------------------------------------------

    assert(msg->size() == m_K);

    // place msg into information positions
    vector<bool> u(m_N, 0);
    for (uint32_t i = 0; i < m_K; i++) {
        u[m_I[i]] = (*msg)[i];
    }

    // encode
    return polar_encode(u);
}
vector<bool> POLAR::sc_decoder(vector<double>* llr)
{
    // make F into bitmap
    vector<bool> F_bit_map(m_N, 0);
    for (auto e : m_F)
        F_bit_map[e] = 1;

    // run the decoding algorithm to decode;
    std::vector<bool> u_cap; // decoded bits
    sc_node_operations(llr, F_bit_map, &u_cap);

    // extract message bits
    vector<bool> msg_cap; msg_cap.reserve(m_K);
    for (auto e : m_I)
        msg_cap.push_back(u_cap[e]);

    return msg_cap;
}
vector<bool> POLAR::scl_decoder(vector<double>* llr, vector<bool> crcG, unsigned nL)
{
    //----------------------------------------------------------------------------------------
    // input:
    // 	   llr: channel info as log likelihood of received bits
    // 	   crcG: crc generator
    // 	   nL: number of list decoders
    // return:
    // 	   decoded message bits
    // Dr J Mao Sep 2021
    //----------------------------------------------------------------------------------------

    uint32_t A = m_K - crcG.size() + 1;  // information bit length excluding padded CRC

    // make F into bitmap
    vector<bool> F_bit_map(m_N, 0);
    for (auto e : m_F) F_bit_map[e] = 1;

    //path metric with first element being 0, the rest are infinity
    vector<double> PM(nL, std::numeric_limits<double>::infinity());
    PM[0] = 0;

    // llr list initialization
    vector<vector<double>> LLR(nL, *llr);

    // hard decisions and coded bits for corresponding llrs
    vector<vector<bool>> hard_decision(nL), u_cap_list(nL);

    // SCL decoding
    scl_node_operations(&LLR, &hard_decision, F_bit_map, &PM, &u_cap_list);

    // CRC checksum to find the best decoding result in the list
    vector<bool> msg_cap(m_K), tmp_cap(m_K);
    vector<size_t> sorted_idx = sort_indexes(PM);

    uint8_t selIdx = sorted_idx[0];
    for (unsigned j = 0; j < nL; j++) {
        for (unsigned i = 0; i < m_I.size(); i++) {
            tmp_cap[i] = u_cap_list[sorted_idx[j]][m_I[i]];
        }
        if (crc_check_sum(&tmp_cap, crcG)) {
            selIdx = sorted_idx[j];
            break;
        }
    }

    for (unsigned i = 0; i < A; i++) {
        msg_cap[i] = u_cap_list[selIdx][m_I[i]];
    }

    return msg_cap;
}
vector<bool> POLAR::rate_matching(vector<bool>* in)
{
    // note, P must be a sorted sequence
    unsigned iL = in->size();
    unsigned pL = m_P.size();

    vector<bool> out(iL - pL);

    unsigned i = 0, j = 0, k = 0; // index to traverse in, P and out

    while (k < iL - pL && j < pL) {
        if (i == m_P[j]) {
            j++; i++;
        }
        else {
            out[k++] = (*in)[i++];
        }
    }

    while (k < iL - pL) {
        out[k++] = (*in)[i++];
    }

    return out;
}
vector<double>  POLAR::rate_recovery(vector<double>* in)
{
    unsigned iL = in->size();
    unsigned pL = m_P.size();

    vector<double> out(iL + pL, 0);

    unsigned i = 0, j = 0, k = 0; // index to traverse in, P and out

    while (i < iL && j < pL) {
        if (k == m_P[j]) {
            k++;
            j++;
        }
        else {
            out[k++] = (*in)[i++];
        }
    }

    while (i < iL) {
        out[k++] = (*in)[i++];
    }
    return out;
}

//Private Core functions-----------------------------------------------------------------------------------------------------------------
vector<bool> POLAR::polar_encode(vector<bool> u)
{
    //---------------------------------------------------------------------------------
    // This function encodes u into x. x = u * F_N
    // where N = 2 ^ n bits,  F_N = F kroncker N, F = [1 0; 1 1];
    // Dr. J Mao 2021 Sep
    //---------------------------------------------------------------------------------

    uint32_t N = u.size();
    if (N == 1)
        return u;
    else {
        // butterfly operation
        vector<bool> u1_xor_u2; u1_xor_u2.reserve(N / 2);
        for (int i = 0; i < N / 2; i++)
            u1_xor_u2.push_back(u[i] xor u[i + N / 2]);

        // recursive encoding for left half
        vector<bool> x; x.reserve(N);
        x = polar_encode(u1_xor_u2);

        // for right half
        vector<bool> x2 = polar_encode(vector<bool>(u.end() - N / 2, u.end()));
        x.insert(x.end(), x2.begin(), x2.end());

        return x;
    }
}
vector<bool> POLAR::sc_node_operations(vector<double>* alpha, vector<bool> F, vector<bool>* u_cap)
{
    unsigned N = alpha->size();
    vector<bool> beta; beta.reserve(N); // hard decisions for corresponding soft bits

    //std::cout << "alpha = [";
    //for (auto e : alpha)
    //	std::cout << e << " ";
    //std::cout<<"]"<< std::endl;

    if (N == 1) { // leaf nodes
        if (F[0] == 1) {
            beta.push_back(0); // frozen bits always zero
        }
        else {
            beta.push_back(((*alpha)[0] < 0) ? 1 : 0); // threshold detection
        }

        u_cap->push_back(beta[0]); // save decoded bits
    }
    else { // non-leaf nodes
        vector<double> a(vector<double>(alpha->begin(), alpha->begin() + N / 2));
        vector<double> b(vector<double>(alpha->end() - N / 2, alpha->end()));

        // calculate llr for its left child
        vector<double> alpha_left(N / 2);
        for (unsigned i = 0; i < N / 2; i++) {
            //alpha_left[i] = 2.0 * atanhf(tanhf(a[i] / 2.0) * tanhf(b[i] / 2.0)); // slower
            int sign = (1 - 2 * (a[i] * b[i] < 0));
            alpha_left[i] = (std::min(abs(a[i]), abs(b[i])) * sign); // much faster with ignorable loss
        }

        // recursive decoding, pass llr message to its left child and expect for corresping hard decision
        vector<bool> F_left(vector<bool>(F.begin(), F.begin() + N / 2));
        vector<bool> beta_left = sc_node_operations(&alpha_left, F_left, u_cap);

        // calculate llr for its right child
        vector<double> alpha_right(N / 2);
        for (unsigned i = 0; i < N / 2; i++) {
            alpha_right[i] = (beta_left[i] == 0) ? b[i] + a[i] : b[i] - a[i];
        }

        // recursive decoding, pass llr message to its right child and expect for corresping hard decision
        vector<bool> F_right(vector<bool>(F.end() - N / 2, F.end()));
        vector<bool> beta_right = sc_node_operations(&alpha_right, F_right, u_cap);

        // combine hard decision from its both children and return to its parent

        for (int i = 0; i < N / 2; i++) {
            beta.push_back(beta_left[i] xor beta_right[i]);
        }
        beta.insert(beta.end(), beta_right.begin(), beta_right.end());
    }
    //std::cout << "beta = [";
    //for (auto e : beta)
    //	std::cout << e << " ";
    //std::cout << "]" << std::endl;
    return beta;
}
void POLAR::scl_node_operations(vector<vector<double>>* llr, vector<vector<bool>>* hard_decision, vector<bool> F, vector<double>* PM, vector<vector<bool>>* u_cap)
{
    //----------------------------------------------------------------------------------------
    // input:
    // 	   llr: channel info as log likelihood of received bits
    // 	   hard decisions: the correspondingly hard bits of the llr, will be caculted as messgaes to pass around nodes
    // 	   F: frozen bits in bitmask format
    // 	   PM: path metric
    // 	   u_cap: decoded bits
    // Dr J Mao Sep 2021
    //----------------------------------------------------------------------------------------

    // the number of bits that the current node carries
    unsigned nb = F.size();
    unsigned nL = llr->size(); // number of decoders;

    // for leaf node
    if (nb == 1) {
        // retrieve the last element of llr list as the soft bit(alpha) of the node
        vector<double> alpha;
        for (unsigned i = 0; i < nL; i++) {
            alpha.push_back((*llr)[i].back());
        }

        // leaf node is frozen bit, hard decision is 0
        if (F[0] == 1) {
            for (unsigned i = 0; i < nL; i++) {
                (*u_cap)[i].push_back(0);
                (*hard_decision)[i].push_back(0);
                // PM update if llr < 0 , no sort and pruning in this case
                (*PM)[i] = (*PM)[i] + abs(alpha[i]) * double((alpha[i] < 0));
            }
        }
        else { // leaf node is information bit
            vector<bool> candi_dec; // candidate decoded bits
            vector<double> candi_PM = *PM; // candidate path metrics
            for (unsigned i = 0; i < nL; i++) { // hard decisions based on llr
                candi_dec.push_back(alpha[i] < 0);
            }
            for (unsigned i = 0; i < nL; i++) { // invert hard decisions and put penalty based on llr
                candi_dec.push_back(alpha[i] >= 0);
                candi_PM.push_back((*PM)[i] + abs(alpha[i]));
            }

            // sort candidates based on path metrics
            vector<size_t> sorted_idx = sort_indexes(candi_PM);
            vector<uint8_t> selected_idx; selected_idx.reserve(nL);

            // sort and prune decoders
            selected_idx.insert(selected_idx.end(), sorted_idx.begin(), sorted_idx.begin() + nL);
            for (auto& ele : selected_idx) {
                if (ele >= nL) ele -= nL;
            }
            sort_and_prune(llr, selected_idx);
            sort_and_prune(hard_decision, selected_idx);
            sort_and_prune(u_cap, selected_idx);

            // update PM and add the newly decoded bits
            for (unsigned i = 0; i < nL; i++) {
                unsigned idx = sorted_idx[i];
                (*PM)[i] = candi_PM[idx];
                (*u_cap)[i].push_back(candi_dec[idx]);
                (*hard_decision)[i].push_back(candi_dec[idx]);
            }
        }
    }
    else {
        // extract the last  node_size bits as the node's sot bits (alpha)
        vector<vector<double>> alpha(nL);
        for (unsigned i = 0; i < nL; i++) {
            alpha[i] = vector<double>((*llr)[i].end() - nb, (*llr)[i].end());
        }

        //check node operation f(a,b) = sgn(a) * sgn(b) * min(|a|,|b|) (min sum approx)
        for (unsigned i = 0; i < nL; i++) {
            for (unsigned j = 0; j < nb / 2; j++) {
                double a = alpha[i][j];
                double b = alpha[i][j + nb / 2];
                int sign = 1 - 2 * (a * b < 0);
                (*llr)[i].push_back(std::min(abs(a), abs(b)) * sign);
            }
        }

        vector<bool> F_left(vector<bool>(F.begin(), F.begin() + nb / 2));
        scl_node_operations(llr, hard_decision, F_left, PM, u_cap);

        //hard decisions returned by left child
        vector<vector<bool>> beta_left(nL);
        for (unsigned i = 0; i < nL; i++) {
            beta_left[i] = vector<bool>((*hard_decision)[i].end() - nb / 2, (*hard_decision)[i].end());
        }

        //extract alpha again as it may have been sorted by the left child
        for (unsigned i = 0; i < nL; i++) {
            alpha[i] = vector <double>((*llr)[i].end() - nb, (*llr)[i].end());
        }

        //bit node handling
        for (unsigned i = 0; i < nL; i++) {
            for (unsigned j = 0; j < nb / 2; j++) {
                double a = alpha[i][j];
                double b = alpha[i][j + nb / 2];
                double c = beta_left[i][j];

                (*llr)[i].push_back(b + (1 - 2 * c) * a);
            }
        }

        // right child decoding
        vector<bool> F_right(vector<bool>(F.end() - nb / 2, F.end()));
        scl_node_operations(llr, hard_decision, F_right, PM, u_cap);

        // butterfly processing does beta = [beta_left xor beta_right, beta_right]
        int start_pos = (*hard_decision)[0].size() - nb;
        for (auto& E : *hard_decision) {
            for (unsigned j = start_pos; j < start_pos + nb / 2; j++) {
                E[j] = (E[j] + E[j + nb / 2]) % 2;
            }
        }
    }

    // delete soft bits from the 2-D llr vectors (equivalent to poping out the calling stack)
    for (auto& ELE : *llr) {
        ELE.erase(ELE.end() - nb, ELE.end());
    }
}
template <typename T>
inline void POLAR::sort_and_prune(vector<vector<T>>* V, vector<uint8_t> P)
{
    vector<vector<T>> tmp;

    for (auto ele : P) {
        tmp.push_back((*V)[ele]);
    }

    V->swap(tmp);
}

// static functions-----------------------------------------------------------------------------------------------------------------------
vector<double> POLAR::channel_polarization_huawei_approx(unsigned N)
{
    //--------------------------------------------------------------------------------
    // return m_N subchannel's weight as a measure of their reliability
    // [ref] 3GPP R1 - 167209 Polar code design and rate matching
    //---------------------------------------------------------------------------------
    vector<double> W(N, 0);

    unsigned n = static_cast<unsigned>(log2(N));

    for (unsigned j = 0; j < N; j++) {
        unsigned tmp = j;

        for (unsigned k = 0; k < n; k++)
        {
            double bit_val = tmp % 2;
            W[j] += bit_val * pow(2, (k * 0.25));
            tmp = tmp / 2;
            if (tmp == 0)
                break;
        }
    }

    return W;
}
template<typename T>
inline vector<size_t> POLAR::sort_indexes(const vector<T>& v)
{
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values
    stable_sort(idx.begin(), idx.end(),
                [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

    return idx;
}

// crc related
bool POLAR::crc_check_sum(vector<bool>* msg, vector<bool> crc_g)
{
    unsigned crc_len = crc_g.size() - 1;
    uint32_t n = msg->size() - crc_len;

    vector<int> msg_crc(msg->begin(), msg->end());
    vector<int> gen(crc_g.begin(), crc_g.end());

    crc_division(&msg_crc, gen);

    for (unsigned i = 0; i < crc_len; i++) {
        if (msg_crc[n + i])
            return 0;
    }

    return 1;
}
void POLAR::crc_division(vector<int>* pad_msg, vector<int> gen)
{
    unsigned r = gen.size() - 1;
    unsigned n = (*pad_msg).size() - r;

    for (unsigned i = 0; i < n; i++)
    {
        if (gen[0] == (*pad_msg)[i])
        {
            for (unsigned j = 0, k = i; j < r + 1; j++, k++) {
                if (!((*pad_msg)[k] ^ gen[j])) {
                    (*pad_msg)[k] = 0;
                }
                else {
                    (*pad_msg)[k] = 1;
                }
            }
        }
    }
}
vector<bool> POLAR::crc_gen(vector<bool>* msg, vector<bool> crc_g)
{
    unsigned crc_len = crc_g.size() - 1;

    vector<int> pad_msg(msg->begin(), msg->end());
    for (uint8_t i = 0; i < crc_len; i++) {
        pad_msg.push_back(0);
    }
    vector<int> gen(crc_g.begin(), crc_g.end());

    crc_division(&pad_msg, gen);

    vector<bool>  crc(pad_msg.end() - crc_len, pad_msg.end());
    return crc;
}
vector<bool> POLAR::crc_generator(const string crc_type)
{
    if (!crc_type.compare("24A"))
        return { 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1 };
    else if (!crc_type.compare("24B"))
        return{ 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1 };
    else if (!crc_type.compare("24C"))
        return{ 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1 };
    else if (!crc_type.compare("16"))
        return{ 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    else if (!crc_type.compare("11"))
        return{ 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    else if (!crc_type.compare("6"))
        return{ 1, 1, 0, 0, 0, 0, 1 };
    else if (!crc_type.compare("1"))
        return{ 1 };
    else
        return {1};
}