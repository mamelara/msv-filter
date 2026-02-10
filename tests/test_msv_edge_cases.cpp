/*******************************************************************************
 * File: tests/test_msv_edge_cases.cpp
 * Description: Edge case tests for MSV filter implementation.
 * Tests boundary conditions, minimal sequences/models, and special cases.
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
float compute_msv(const DigitalResidue* digital_sequence, int sequence_length,
                  const HMMProfile& profile, DPMatrix& dp_matrix, float expected_hit_count);

namespace msv_test {

// ============================================================================
// Edge Case Test Fixtures
// ============================================================================

// Minimal model (M=1) and minimal sequence (L=1)
struct MinimalTestCase {
    static constexpr int MODEL_LENGTH = 1;
    static constexpr int SEQUENCE_LENGTH = 1;
    static constexpr float MATCH_SCORE = 1.0f;
    static constexpr float EXPECTED_SCORE = 1.0f;
    
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// Empty sequence (L=0) - should handle gracefully
struct EmptySequenceTestCase {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 0;
    static constexpr float MATCH_SCORE = 1.0f;
    static constexpr float EXPECTED_SCORE = 0.0f;  // No sequence = no score
    
    static std::vector<DigitalResidue> get_sequence() {
        // Just sentinels, no actual sequence data
        std::vector<DigitalResidue> seq(2);
        seq[0] = digitalResidueSentinel;
        seq[1] = digitalResidueSentinel;
        return seq;
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// Empty model (M=0) - should handle gracefully
struct EmptyModelTestCase {
    static constexpr int MODEL_LENGTH = 0;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float MATCH_SCORE = 1.0f;
    static constexpr float EXPECTED_SCORE = 0.0f;  // No model = no score
    
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D, RES_E, RES_F});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        // Create profile with M=1 but mark actual model_length as 0
        HMMProfile profile(1, &abc);
        profile.model_length = 0;
        profile.name = "empty_model";
        profile.max_length = 0;
        return profile;
    }
    
    static DPMatrix get_dp_matrix() {
        // DP matrix with M=0 is unusual, use minimal allocation
        return DPMatrix(0, SEQUENCE_LENGTH);
    }
};

// Very large negative scores - tests numerical stability
struct LargeNegativeScoresTestCase {
    static constexpr int MODEL_LENGTH = 3;
    static constexpr int SEQUENCE_LENGTH = 3;
    static constexpr float GOOD_SCORE = 5.0f;
    static constexpr float BAD_SCORE = -100.0f;
    static constexpr float EXPECTED_SCORE = 5.0f;  // Only the good position
    
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_A, RES_A});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        HMMProfile profile(MODEL_LENGTH, &abc);
        profile.model_length = MODEL_LENGTH;
        profile.name = "large_negative_model";
        profile.max_length = 100;
        
        // Position 1: A=5.0, others=-100
        // Position 2: all=-100
        // Position 3: all=-100
        for (int k = 1; k <= MODEL_LENGTH; k++) {
            for (int x = 0; x < abc.K; x++) {
                if (k == 1 && x == RES_A) {
                    profile.match_score(k, x) = GOOD_SCORE;
                } else {
                    profile.match_score(k, x) = BAD_SCORE;
                }
            }
        }
        
        return profile;
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// Very large positive scores - tests numerical stability
struct LargePositiveScoresTestCase {
    static constexpr int MODEL_LENGTH = 3;
    static constexpr int SEQUENCE_LENGTH = 3;
    static constexpr float SCORE = 1000.0f;
    static constexpr float EXPECTED_SCORE = 3000.0f;  // 3 * 1000
    
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// All negative scores - MSV should find best (least negative) segment
struct AllNegativeScoresTestCase {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float SCORE = -2.0f;
    // With all negative, best is to take single best position: -2.0
    // or possibly 0.0 if empty alignment is allowed
    static constexpr float EXPECTED_SCORE = -2.0f;
    
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D, RES_E, RES_F});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// Zero scores everywhere
struct ZeroScoresTestCase {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float SCORE = 0.0f;
    static constexpr float EXPECTED_SCORE = 0.0f;
    
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D, RES_E, RES_F});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// Long model, short sequence (boundary test)
struct LongModelShortSequenceTestCase {
    static constexpr int MODEL_LENGTH = 100;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float SCORE = 2.0f;
    // Can only align 5 positions due to sequence limit
    static constexpr float EXPECTED_SCORE = 10.0f;  // 5 * 2.0
    
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D, RES_E, RES_F});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// Short model, long sequence (boundary test)
struct ShortModelLongSequenceTestCase {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 100;
    static constexpr float SCORE = 2.0f;
    // Can only align 5 positions due to model limit
    static constexpr float EXPECTED_SCORE = 10.0f;  // 5 * 2.0
    
    static std::vector<DigitalResidue> get_sequence() {
        std::vector<DigitalResidue> residues;
        residues.reserve(100);
        for (int i = 0; i < 100; i++) {
            residues.push_back(static_cast<DigitalResidue>(i % 20));
        }
        return create_digital_sequence(residues);
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// Sequence with degenerate/ambiguous residues
struct DegenerateResiduesTestCase {
    static constexpr int MODEL_LENGTH = 3;
    static constexpr int SEQUENCE_LENGTH = 3;
    static constexpr float EXPECTED_SCORE = 3.0f;
    
    static std::vector<DigitalResidue> get_sequence() {
        // Use standard residues for this test
        return create_digital_sequence({RES_A, RES_C, RES_D});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        HMMProfile profile(MODEL_LENGTH, &abc);
        profile.model_length = MODEL_LENGTH;
        profile.name = "degenerate_test_model";
        profile.max_length = 100;
        
        // Simple scoring: all positions score 1.0 for all residues
        for (int k = 1; k <= MODEL_LENGTH; k++) {
            for (int x = 0; x < abc.K; x++) {
                profile.match_score(k, x) = 1.0f;
            }
        }
        
        return profile;
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

} // namespace msv_test

// ============================================================================
// Test Fixture for Edge Case Tests
// ============================================================================
class MSVEdgeCaseTest : public ::testing::Test {
protected:
    static const AminoAcidAlphabet* alphabet;
    
    static void SetUpTestSuite() {
        alphabet = &msv_test::get_test_alphabet();
    }
    
    template<typename TestCase>
    void run_edge_case_test() {
        std::vector<DigitalResidue> digital_sequence = TestCase::get_sequence();
        HMMProfile profile = TestCase::get_profile(*alphabet);
        DPMatrix dp_matrix = TestCase::get_dp_matrix();
        
        float actual_score = compute_msv(
            digital_sequence.data(),
            TestCase::SEQUENCE_LENGTH,
            profile,
            dp_matrix,
            1.0f
        );
        
        ASSERT_NEAR(TestCase::EXPECTED_SCORE, actual_score, 0.001f)
            << "MSV score mismatch for edge case: " << profile.name;
    }
};

const AminoAcidAlphabet* MSVEdgeCaseTest::alphabet = nullptr;

// ============================================================================
// Edge Case Tests
// ============================================================================

// Minimal model and sequence (M=1, L=1)
// This tests the smallest non-trivial case
TEST_F(MSVEdgeCaseTest, MinimalCase) {
    run_edge_case_test<msv_test::MinimalTestCase>();
}

// Empty sequence (L=0)
// The function should return 0.0 when there's no sequence to align
TEST_F(MSVEdgeCaseTest, EmptySequence) {
    run_edge_case_test<msv_test::EmptySequenceTestCase>();
}

// Empty model (M=0)
// The function should return 0.0 when there's no model to align to
TEST_F(MSVEdgeCaseTest, EmptyModel) {
    run_edge_case_test<msv_test::EmptyModelTestCase>();
}

// Large negative scores mixed with positive
// Tests that the algorithm correctly handles extreme score differences
TEST_F(MSVEdgeCaseTest, LargeNegativeScores) {
    run_edge_case_test<msv_test::LargeNegativeScoresTestCase>();
}

// Large positive scores
// Tests numerical stability with large values
TEST_F(MSVEdgeCaseTest, LargePositiveScores) {
    run_edge_case_test<msv_test::LargePositiveScoresTestCase>();
}

// All negative scores
// MSV should find the best (least negative) single position
TEST_F(MSVEdgeCaseTest, AllNegativeScores) {
    run_edge_case_test<msv_test::AllNegativeScoresTestCase>();
}

// All zero scores
// Every alignment should give 0.0
TEST_F(MSVEdgeCaseTest, ZeroScores) {
    run_edge_case_test<msv_test::ZeroScoresTestCase>();
}

// Long model (M=100) with short sequence (L=5)
// Tests that the sequence length is the limiting factor
TEST_F(MSVEdgeCaseTest, LongModelShortSequence) {
    run_edge_case_test<msv_test::LongModelShortSequenceTestCase>();
}

// Short model (M=5) with long sequence (L=100)
// Tests that the model length is the limiting factor
TEST_F(MSVEdgeCaseTest, ShortModelLongSequence) {
    run_edge_case_test<msv_test::ShortModelLongSequenceTestCase>();
}

// Degenerate/ambiguous residue handling
TEST_F(MSVEdgeCaseTest, DegenerateResidues) {
    run_edge_case_test<msv_test::DegenerateResiduesTestCase>();
}

// ============================================================================
// Sentinel Verification Tests
// ============================================================================
// These tests verify that sentinels are properly handled

TEST_F(MSVEdgeCaseTest, VerifySentinelsUnchanged) {
    using namespace msv_test;
    
    auto digital_sequence = MinimalTestCase::get_sequence();
    
    // Store original sentinel values
    DigitalResidue sentinel_before_0 = digital_sequence[0];
    DigitalResidue sentinel_before_end = digital_sequence[digital_sequence.size() - 1];
    
    HMMProfile profile = MinimalTestCase::get_profile(*alphabet);
    DPMatrix dp_matrix = MinimalTestCase::get_dp_matrix();
    
    // Call compute_msv
    compute_msv(digital_sequence.data(), 1, profile, dp_matrix, 1.0f);
    
    // Verify sentinels are unchanged
    EXPECT_EQ(sentinel_before_0, digital_sequence[0]);
    EXPECT_EQ(sentinel_before_end, digital_sequence[digital_sequence.size() - 1]);
    EXPECT_EQ(digitalResidueSentinel, digital_sequence[0]);
    EXPECT_EQ(digitalResidueSentinel, digital_sequence[digital_sequence.size() - 1]);
}

// ============================================================================
// Boundary Value Tests
// ============================================================================

TEST_F(MSVEdgeCaseTest, SingleResidueModelMultiplePositions) {
    // M=1, L=10 - single position model, longer sequence
    // Should only be able to score one position
    std::vector<DigitalResidue> seq = msv_test::create_digital_sequence({
        msv_test::RES_A, msv_test::RES_C, msv_test::RES_D, msv_test::RES_E, msv_test::RES_F,
        msv_test::RES_G, msv_test::RES_H, msv_test::RES_I, msv_test::RES_K, msv_test::RES_L
    });
    
    HMMProfile profile(1, alphabet);
    profile.model_length = 1;
    profile.name = "single_pos_long_seq";
    profile.max_length = 100;
    
    // Position 1 scores 5.0 for all residues
    for (int x = 0; x < alphabet->K; x++) {
        profile.match_score(1, x) = 5.0f;
    }
    
    DPMatrix dp_matrix(1, 10);
    
    float score = compute_msv(seq.data(), 10, profile, dp_matrix, 1.0f);
    
    // Can only align 1 position, score = 5.0
    EXPECT_NEAR(5.0f, score, 0.001f);
}

TEST_F(MSVEdgeCaseTest, ModelAndSequenceBothSingle) {
    // M=1, L=1 - absolute minimum
    std::vector<DigitalResidue> seq = msv_test::create_digital_sequence({msv_test::RES_A});
    
    HMMProfile profile(1, alphabet);
    profile.model_length = 1;
    profile.name = "minimal";
    profile.max_length = 100;
    profile.match_score(1, msv_test::RES_A) = 7.5f;
    
    DPMatrix dp_matrix(1, 1);
    
    float score = compute_msv(seq.data(), 1, profile, dp_matrix, 1.0f);
    
    EXPECT_NEAR(7.5f, score, 0.001f);
}

// ============================================================================
// Memory Boundary Tests
// ============================================================================
// These tests verify that the DP matrix boundaries are respected

TEST_F(MSVEdgeCaseTest, DPMatrixDimensions) {
    using namespace msv_test;
    
    auto seq = ConstantAllOnesTest::get_sequence();
    auto profile = ConstantAllOnesTest::get_profile(*alphabet);
    auto dp_matrix = ConstantAllOnesTest::get_dp_matrix();
    
    // Verify DP matrix has correct dimensions
    EXPECT_EQ(5, dp_matrix.model_length);
    EXPECT_EQ(5, dp_matrix.sequence_length);
    EXPECT_EQ(6, dp_matrix.dp.size());  // L+1 rows (0..5)
    
    // Each row should have (M+1) * 3 cells
    // M+1 = 6 positions (0..5), 3 states each = 18 cells per row
    EXPECT_EQ(18, dp_matrix.dp[0].size());
}
