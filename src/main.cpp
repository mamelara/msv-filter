/*******************************************************************************
 * File: src/main.cpp
 * Description: Demonstrates mocked inputs for p7_GMSV function.
 * 
 * This file creates all the necessary inputs to call a re-implemented
 * p7_GMSV function with mocked data:
 *   - DigitalResidue digital sequence
 *   - P7_PROFILE structure
 *   - P7_GMX DP matrix
 *   - expected_hit_count parameter
 *   - msv_score return value
 ******************************************************************************/

#include <iostream>
#include <vector>
#include <cmath>
#include "hmmer_types.hpp"
#include "aa_alphabet.hpp"
#include "profile.hpp"
#include "dp_matrix.hpp"
#include "mock_data.hpp"

/*******************************************************************************
 * Example signature of the MSV function to be implemented:
 * 
 * int p7_GMSV(const DigitalResidue *digital_sequence, int sequence_length, const P7_PROFILE *gm, 
 *             P7_GMX *gx, float expected_hit_count, float *msv_score)
 * 
 * Where:
 *   - digital_sequence: Digital sequence (1-indexed array of residue indices)
 *   - sequence_length:   Sequence length
 *   - gm:  Profile with match scores (MSV only uses MSC, not transitions)
 *   - gx:  DP matrix with match states and special states
 *   - expected_hit_count:  Expected number of hits (typically 2.0)
 *   - msv_score: Return value for MSV score
 ******************************************************************************/

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "MSV Filter - Mock Input Generator" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // --- Step 1: Create Alphabet ---
    std::cout << "\n[1] Creating Amino Acid Alphabet..." << std::endl;
    AminoAcidAlphabet abc;
    std::cout << "    Alphabet size (K): " << abc.K << std::endl;
    std::cout << "    Total symbols (Kp): " << abc.Kp << std::endl;
    std::cout << "    Symbols: " << abc.sym << std::endl;
    
    // --- Step 2: Create Mock Digital Sequence ---
    std::cout << "\n[2] Creating Mock Digital Sequence..." << std::endl;
    int sequence_length = 15;  // Sequence length
    std::vector<DigitalResidue> digital_sequence = MockDataGenerator::create_simple_sequence(sequence_length, abc);
    std::cout << "    Length (sequence_length): " << sequence_length << std::endl;
    MockDataGenerator::print_sequence(digital_sequence, sequence_length, abc);
    std::cout << "    Digital representation (first 10): ";
    for (int i = 1; i <= std::min(10, sequence_length); i++) {
        std::cout << (int)digital_sequence[i] << " ";
    }
    std::cout << std::endl;
    
    // --- Step 3: Create Mock Profile ---
    std::cout << "\n[3] Creating Mock Profile..." << std::endl;
    int model_length = 10;  // Model length
    HMMProfile profile = MockDataGenerator::create_simple_profile(model_length, abc);
    std::cout << "    Model length (model_length): " << profile.model_length << std::endl;
    std::cout << "    Max length: " << profile.max_length << std::endl;
    std::cout << "    Match scores (first 3 positions, first 5 residues):" << std::endl;
    for (int k = 1; k <= std::min(3, model_length); k++) {
        std::cout << "      k=" << k << ": ";
        for (int x = 0; x < std::min(5, abc.K); x++) {
            std::cout << profile.match_score(k, x) << " ";
        }
        std::cout << std::endl;
    }
    
    // --- Step 4: Create DP Matrix ---
    std::cout << "\n[4] Creating DP Matrix..." << std::endl;
    DPMatrix dp_matrix(model_length, sequence_length);
    std::cout << "    Matrix dimensions: model_length=" << dp_matrix.model_length << ", sequence_length=" << dp_matrix.sequence_length << std::endl;
    std::cout << "    DP cells: " << (sequence_length + 1) << " rows x " << ((model_length + 1) * 3) << " cols" << std::endl;
    std::cout << "    Special states: " << (sequence_length + 1) << " rows x " << 5 << " cols" << std::endl;
    
    // --- Step 5: Set up other parameters ---
    std::cout << "\n[5] Setting up additional parameters..." << std::endl;
    float expected_hit_count = 2.0f;  // Expected number of hits
    float msv_score = 0.0f;  // Will hold the result
    std::cout << "    expected_hit_count (expected hits): " << expected_hit_count << std::endl;
    std::cout << "    msv_score (output): " << msv_score << " (will be set by MSV)" << std::endl;
    
    // --- Step 6: Show memory layout ---
    std::cout << "\n[6] Memory Layout for p7_GMSV inputs:" << std::endl;
    std::cout << "    digital_sequence (DigitalResidue*): " << digital_sequence.size() << " bytes" << std::endl;
    std::cout << "      - Index 0: " << (int)digital_sequence[0] << " (sentinel)" << std::endl;
    std::cout << "      - Index 1.." << sequence_length << ": residues" << std::endl;
    std::cout << "      - Index " << (sequence_length+1) << ": " << (int)digital_sequence[sequence_length+1] << " (sentinel)" << std::endl;
    
    std::cout << "\n    gm (P7_PROFILE*): " << std::endl;
    std::cout << "      - model_length: " << profile.model_length << std::endl;
    std::cout << "      - tsc: " << profile.tsc.size() << " floats (transitions)" << std::endl;
    std::cout << "      - rsc: " << profile.rsc.size() << " x " << profile.rsc[0].size() << " floats (emissions)" << std::endl;
    std::cout << "      - xsc: " << p7P_NXSTATES << " x " << p7P_NXTRANS << " floats (special transitions)" << std::endl;
    
    std::cout << "\n    gx (P7_GMX*): " << std::endl;
    std::cout << "      - model_length: " << dp_matrix.model_length << ", sequence_length: " << dp_matrix.sequence_length << std::endl;
    std::cout << "      - dp: " << dp_matrix.dp.size() << " rows x " << dp_matrix.dp[0].size() << " cols" << std::endl;
    std::cout << "      - xmx: " << dp_matrix.xmx.size() << " cells" << std::endl;
    
    // --- Step 7: Summary ---
    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary: Ready to call p7_GMSV" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Function signature:" << std::endl;
    std::cout << "  int p7_GMSV(const DigitalResidue *digital_sequence, int sequence_length, const P7_PROFILE *gm," << std::endl;
    std::cout << "              P7_GMX *gx, float expected_hit_count, float *msv_score)" << std::endl;
    std::cout << std::endl;
    std::cout << "Arguments ready:" << std::endl;
    std::cout << "  - digital_sequence: &digital_sequence[0] (1-indexed, 0 is sentinel)" << std::endl;
    std::cout << "  - sequence_length: " << sequence_length << std::endl;
    std::cout << "  - gm: &profile" << std::endl;
    std::cout << "  - gx: &dp_matrix" << std::endl;
    std::cout << "  - expected_hit_count: " << expected_hit_count << std::endl;
    std::cout << "  - msv_score: &msv_score" << std::endl;
    
    std::cout << "\nNote: MSV algorithm only uses:" << std::endl;
    std::cout << "  - gm->rsc[residue][k * 2 + 0] (match scores)" << std::endl;
    std::cout << "  - gm->model_length (model length)" << std::endl;
    std::cout << "  - gx->dp[i][k * 3 + 0] (match states)" << std::endl;
    std::cout << "  - gx->xmx[i * 5 + s] (special states: E,N,J,B,C)" << std::endl;
    
    return 0;
}
