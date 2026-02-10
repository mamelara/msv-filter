/*******************************************************************************
 * File: tests/test_msv_basic.cpp
 * Description: Basic functionality tests for MSV filter implementation.
 * Uses hardcoded test vectors from test_vectors.hpp.
 ******************************************************************************/

#include <gtest/gtest.h>
#include "test_vectors.hpp"
#include "hmmer_types.hpp"
#include "profile.hpp"
#include "dp_matrix.hpp"
#include "mock_data.hpp"
#include "aa_alphabet.hpp"

// ============================================================================
// Forward Declaration - User will implement this function
// ============================================================================
// The MSV (Multiple Segment Viterbi) filter computes the maximum scoring
// segment alignment between a sequence and an HMM profile.
// 
// Parameters:
//   digital_sequence: 1-indexed digital sequence with sentinels at 0 and L+1
//   sequence_length:  Length of the actual sequence (sentinels excluded)
//   profile:          HMM profile with match scores
//   dp_matrix:        Dynamic programming matrix for storing intermediate values
//   expected_hit_count: Expected number of hits (for normalization)
//
// Returns:
//   The MSV score (maximum scoring segment score)
//
// Note: This is a placeholder - the user will implement the actual algorithm
float compute_msv(const DigitalResidue* digital_sequence, int sequence_length,
                  const HMMProfile& profile, DPMatrix& dp_matrix, float expected_hit_count);

// ============================================================================
// Test Fixture for MSV Basic Tests
// ============================================================================
class MSVBasicTest : public ::testing::Test {
protected:
    // Alphabet is initialized once for all tests
    static const AminoAcidAlphabet* alphabet;
    
    static void SetUpTestSuite() {
        alphabet = &msv_test::get_test_alphabet();
    }
    
    // Helper to run a test case and verify result
    template<typename TestCase>
    void run_test_case() {
        // Get test data from test case
        std::vector<DigitalResidue> digital_sequence = TestCase::get_sequence();
        HMMProfile profile = TestCase::get_profile(*alphabet);
        DPMatrix dp_matrix = TestCase::get_dp_matrix();
        
        // Call the MSV function (placeholder)
        float actual_score = compute_msv(
            digital_sequence.data(), 
            TestCase::SEQUENCE_LENGTH,
            profile, 
            dp_matrix, 
            1.0f  // expected_hit_count
        );
        
        // Verify with tolerance for floating point
        ASSERT_NEAR(TestCase::EXPECTED_SCORE, actual_score, 0.001f)
            << "MSV score mismatch for test case: " << profile.name;
    }
};

const AminoAcidAlphabet* MSVBasicTest::alphabet = nullptr;

// ============================================================================
// Test Case 1: ConstantAllOnes
// ============================================================================
// All match scores = 1.0, M=5, L=5
// Expected calculation:
//   With all scores equal to 1.0, MSV should find the best segment
//   Since the model and sequence have the same length (5), 
//   the optimal alignment matches all 5 positions.
//   Score = 5 positions * 1.0 = 5.0
// ============================================================================
TEST_F(MSVBasicTest, ConstantAllOnes) {
    run_test_case<msv_test::ConstantAllOnesTest>();
}

// ============================================================================
// Test Case 2: ConstantAllTwos
// ============================================================================
// All match scores = 2.0, M=5, L=5
// Expected calculation:
//   Similar to ConstantAllOnes, but with doubled scores.
//   Score = 5 positions * 2.0 = 10.0
// ============================================================================
TEST_F(MSVBasicTest, ConstantAllTwos) {
    run_test_case<msv_test::ConstantAllTwosTest>();
}

// ============================================================================
// Test Case 3: SinglePositionModel
// ============================================================================
// M=1 (single position model), L=5
// Expected calculation:
//   The model has only one position, so MSV can only match one residue
//   regardless of sequence length.
//   Score = 1 position * 1.0 = 1.0
// ============================================================================
TEST_F(MSVBasicTest, SinglePositionModel) {
    run_test_case<msv_test::SinglePositionModelTest>();
}

// ============================================================================
// Test Case 4: SingleResidueSequence
// ============================================================================
// M=5, L=1 (single residue sequence)
// Expected calculation:
//   The sequence has only one residue, so MSV can match at most one position.
//   Score = 1 residue * 1.0 = 1.0
// ============================================================================
TEST_F(MSVBasicTest, SingleResidueSequence) {
    run_test_case<msv_test::SingleResidueSequenceTest>();
}

// ============================================================================
// Test Case 5: AlternatingPattern
// ============================================================================
// M=10, L=10 with position k preferring residue (k % K)
// Match score = 3.0 for preferred, -1.0 for others
// Expected calculation:
//   Sequence is designed to match: residues 0,1,2,3,4,5,6,7,8,9
//   Model position k prefers residue (k-1), so:
//     Position 1 matches residue 0: score = 3.0
//     Position 2 matches residue 1: score = 3.0
//     ...and so on
//   Total score = 10 positions * 3.0 = 30.0
// ============================================================================
TEST_F(MSVBasicTest, AlternatingPattern) {
    run_test_case<msv_test::AlternatingPatternTest>();
}

// ============================================================================
// Test Case 6: AllSameResidue
// ============================================================================
// Sequence "AAAAA" (all Alanine), M=5, L=5
// Only position 1 likes Alanine (score 3.0), others give -1.0
// Expected calculation:
//   MSV finds the maximum scoring segment.
//   Aligning position 1 with any A gives score 3.0.
//   Aligning position 1 and another position gives 3.0 + (-1.0) = 2.0 (worse)
//   Best segment is just position 1: score = 3.0
// ============================================================================
TEST_F(MSVBasicTest, AllSameResidue) {
    run_test_case<msv_test::AllSameResidueTest>();
}

// ============================================================================
// Test Case 7: AllDifferentResidues
// ============================================================================
// Full alphabet sequence: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
// M=20, L=20
// Each position k matches its corresponding residue with score 2.0
// Expected calculation:
//   Perfect match: each position aligned with its preferred residue
//   Score = 20 positions * 2.0 = 40.0
// ============================================================================
TEST_F(MSVBasicTest, AllDifferentResidues) {
    run_test_case<msv_test::AllDifferentResiduesTest>();
}

// ============================================================================
// Test Case 8: ShorterSequence
// ============================================================================
// L=3, M=10 (short sequence, longer model)
// Constant scores = 2.0
// Expected calculation:
//   Sequence is the limiting factor with only 3 residues.
//   MSV can align at most 3 positions.
//   Score = 3 * 2.0 = 6.0
// ============================================================================
TEST_F(MSVBasicTest, ShorterSequence) {
    run_test_case<msv_test::ShorterSequenceTest>();
}

// ============================================================================
// Test Case 9: LongerSequence
// ============================================================================
// L=20, M=5 (long sequence, shorter model)
// Constant scores = 1.5
// Expected calculation:
//   Model is the limiting factor with only 5 positions.
//   MSV can align at most 5 positions.
//   Score = 5 * 1.5 = 7.5
// ============================================================================
TEST_F(MSVBasicTest, LongerSequence) {
    run_test_case<msv_test::LongerSequenceTest>();
}

// ============================================================================
// Test Case 10: MixedScores
// ============================================================================
// M=4, L=4 with specific position-dependent scores
// Position 1: A=2.0, Position 2: C=3.0, Position 3: D=2.0, Position 4: E=3.0
// Sequence: A, C, D, E
// Expected calculation:
//   Position 1 aligns with A: 2.0
//   Position 2 aligns with C: 3.0
//   Position 3 aligns with D: 2.0
//   Position 4 aligns with E: 3.0
//   Total: 2.0 + 3.0 + 2.0 + 3.0 = 10.0
// ============================================================================
TEST_F(MSVBasicTest, MixedScores) {
    run_test_case<msv_test::MixedScoresTest>();
}

// ============================================================================
// Direct Value Verification Tests
// ============================================================================
// These tests verify specific expected values without relying on compute_msv
// Useful for debugging and validating test vector correctness

TEST_F(MSVBasicTest, VerifyTestVectors) {
    // Verify ConstantAllOnes has correct constants
    EXPECT_EQ(5, msv_test::ConstantAllOnesTest::MODEL_LENGTH);
    EXPECT_EQ(5, msv_test::ConstantAllOnesTest::SEQUENCE_LENGTH);
    EXPECT_FLOAT_EQ(1.0f, msv_test::ConstantAllOnesTest::MATCH_SCORE);
    EXPECT_FLOAT_EQ(5.0f, msv_test::ConstantAllOnesTest::EXPECTED_SCORE);
    
    // Verify sequence has sentinels
    auto seq = msv_test::ConstantAllOnesTest::get_sequence();
    EXPECT_EQ(digitalResidueSentinel, seq[0]);
    EXPECT_EQ(digitalResidueSentinel, seq[6]);
    EXPECT_EQ(msv_test::RES_A, seq[1]);
    EXPECT_EQ(msv_test::RES_C, seq[2]);
    
    // Verify profile has correct dimensions
    auto profile = msv_test::ConstantAllOnesTest::get_profile(*alphabet);
    EXPECT_EQ(5, profile.model_length);
    
    // All match scores should be 1.0
    for (int k = 1; k <= 5; k++) {
        for (int x = 0; x < alphabet->K; x++) {
            EXPECT_FLOAT_EQ(1.0f, profile.match_score(k, x))
                << "Mismatch at k=" << k << ", x=" << x;
        }
    }
}

TEST_F(MSVBasicTest, VerifyAlternatingPattern) {
    auto profile = msv_test::AlternatingPatternTest::get_profile(*alphabet);
    
    // Position 1 should prefer residue 0 (A)
    EXPECT_FLOAT_EQ(3.0f, profile.match_score(1, msv_test::RES_A));
    EXPECT_FLOAT_EQ(-1.0f, profile.match_score(1, msv_test::RES_C));
    
    // Position 2 should prefer residue 1 (C)
    EXPECT_FLOAT_EQ(3.0f, profile.match_score(2, msv_test::RES_C));
    EXPECT_FLOAT_EQ(-1.0f, profile.match_score(2, msv_test::RES_A));
}

TEST_F(MSVBasicTest, VerifyAllSameResidue) {
    auto profile = msv_test::AllSameResidueTest::get_profile(*alphabet);
    
    // Position 1 should give 3.0 for Alanine
    EXPECT_FLOAT_EQ(3.0f, profile.match_score(1, msv_test::RES_A));
    
    // Other positions should give -1.0 for Alanine
    EXPECT_FLOAT_EQ(-1.0f, profile.match_score(2, msv_test::RES_A));
    EXPECT_FLOAT_EQ(-1.0f, profile.match_score(3, msv_test::RES_A));
}
