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
#include "hmmer_types.hpp"
#include "aa_alphabet.hpp"
#include "profile.hpp"
#include "dp_matrix.hpp"
#include "mock_data.hpp"
#include "msv_filter.hpp"

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

    AminoAcidAlphabet abc;

    int sequence_length = 15;  // Sequence length
    std::vector<DigitalResidue> digital_sequence = MockDataGenerator::create_simple_sequence(sequence_length, abc);
    MockDataGenerator::print_sequence(digital_sequence, sequence_length, abc);

    // --- Step 3: Create Mock Profile ---
    int model_length = 10;  // Model length
    HMMProfile profile = MockDataGenerator::create_simple_profile(model_length, abc);

    DPMatrix dp_matrix(model_length, sequence_length);

    float expected_hit_count = 2.0f;  // Expected number of hits
    float msv_score = 0.0f;  // Will hold the result

  msv_filter(digital_sequence.data(), sequence_length, profile, dp_matrix, expected_hit_count, &msv_score);

  std::cout << "MSV score: " << msv_score << std::endl;



    return 0;
}
