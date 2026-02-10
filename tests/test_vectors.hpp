/*******************************************************************************
 * File: tests/test_vectors.hpp
 * Description: Hardcoded test vectors for MSV filter unit tests.
 * All test data is deterministic - no random generation.
 * 
 * Test cases include:
 *   - ConstantAllOnes: All match scores = 1.0, M=5, L=5
 *   - ConstantAllTwos: All match scores = 2.0, M=5, L=5
 *   - SinglePositionModel: M=1, L=5
 *   - SingleResidueSequence: M=5, L=1
 *   - AlternatingPattern: Position k prefers residue (k % K)
 *   - AllSameResidue: Sequence "AAAAA"
 *   - AllDifferentResidues: Full alphabet sequence
 *   - ShorterSequence: L=3, M=10
 *   - LongerSequence: L=20, M=5
 ******************************************************************************/

#ifndef MSV_FILTER_TEST_VECTORS_HPP
#define MSV_FILTER_TEST_VECTORS_HPP

#include <vector>
#include <array>
#include <memory>
#include "hmmer_types.hpp"
#include "aa_alphabet.hpp"
#include "profile.hpp"
#include "dp_matrix.hpp"

namespace msv_test {

// ============================================================================
// Test Alphabet Setup
// ============================================================================

// Amino acid alphabet used for all tests
// Standard 20 amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
// Digital encoding: A=0, C=1, D=2, E=3, F=4, G=5, H=6, I=7, K=8, L=9, 
//                   M=10, N=11, P=12, Q=13, R=14, S=15, T=16, V=17, W=18, Y=19

// Helper function to get singleton alphabet instance
inline const AminoAcidAlphabet& get_test_alphabet() {
    static std::unique_ptr<AminoAcidAlphabet> abc;
    if (!abc) {
        abc = std::make_unique<AminoAcidAlphabet>();
    }
    return *abc;
}

// ============================================================================
// Helper Functions for Creating Test Fixtures
// ============================================================================

// Create a digital sequence with sentinels (index 0 and L+1 are sentinels)
// The actual sequence is at indices 1..sequence_length
inline std::vector<DigitalResidue> create_digital_sequence(
    const std::vector<DigitalResidue>& residues) 
{
    int sequence_length = static_cast<int>(residues.size());
    std::vector<DigitalResidue> digital_sequence(sequence_length + 2);
    
    // Set sentinels
    digital_sequence[0] = digitalResidueSentinel;
    digital_sequence[sequence_length + 1] = digitalResidueSentinel;
    
    // Copy residues to positions 1..sequence_length
    for (int i = 0; i < sequence_length; i++) {
        digital_sequence[i + 1] = residues[i];
    }
    
    return digital_sequence;
}

// Create a profile with constant match scores for all residues
inline HMMProfile create_constant_score_profile(
    int model_length, 
    float match_score,
    const AminoAcidAlphabet& abc)
{
    HMMProfile profile(model_length, &abc);
    profile.model_length = model_length;
    profile.name = "constant_score_model";
    profile.max_length = 100;
    
    // Set all match scores to the same value for all positions and residues
    for (int k = 1; k <= model_length; k++) {
        for (int x = 0; x < abc.K; x++) {
            profile.match_score(k, x) = match_score;
        }
    }
    
    return profile;
}

// Create a profile with position-specific preferred residues
// Position k (1-indexed) prefers residue ((k-1) % K)
inline HMMProfile create_alternating_pattern_profile(
    int model_length,
    float match_score,
    float mismatch_score,
    const AminoAcidAlphabet& abc)
{
    HMMProfile profile(model_length, &abc);
    profile.model_length = model_length;
    profile.name = "alternating_pattern_model";
    profile.max_length = 100;
    
    // Each position k prefers residue (k-1) % K
    for (int k = 1; k <= model_length; k++) {
        int preferred_residue = (k - 1) % abc.K;
        for (int x = 0; x < abc.K; x++) {
            if (x == preferred_residue) {
                profile.match_score(k, x) = match_score;    // High score for preferred
            } else {
                profile.match_score(k, x) = mismatch_score; // Lower score for others
            }
        }
    }
    
    return profile;
}

// Create a profile with a specific residue scoring pattern
inline HMMProfile create_specific_pattern_profile(
    int model_length,
    const std::vector<std::vector<float>>& scores_per_position,
    const AminoAcidAlphabet& abc)
{
    HMMProfile profile(model_length, &abc);
    profile.model_length = model_length;
    profile.name = "specific_pattern_model";
    profile.max_length = 100;
    
    for (int k = 1; k <= model_length && k <= static_cast<int>(scores_per_position.size()); k++) {
        for (int x = 0; x < abc.K && x < static_cast<int>(scores_per_position[k-1].size()); x++) {
            profile.match_score(k, x) = scores_per_position[k-1][x];
        }
    }
    
    return profile;
}

// ============================================================================
// Digital Residue Constants (for clarity in test definitions)
// ============================================================================

// Standard amino acid digital encodings
constexpr DigitalResidue RES_A = 0;   // Alanine
constexpr DigitalResidue RES_C = 1;   // Cysteine
constexpr DigitalResidue RES_D = 2;   // Aspartic Acid
constexpr DigitalResidue RES_E = 3;   // Glutamic Acid
constexpr DigitalResidue RES_F = 4;   // Phenylalanine
constexpr DigitalResidue RES_G = 5;   // Glycine
constexpr DigitalResidue RES_H = 6;   // Histidine
constexpr DigitalResidue RES_I = 7;   // Isoleucine
constexpr DigitalResidue RES_K = 8;   // Lysine
constexpr DigitalResidue RES_L = 9;   // Leucine
constexpr DigitalResidue RES_M = 10;  // Methionine
constexpr DigitalResidue RES_N = 11;  // Asparagine
constexpr DigitalResidue RES_P = 12;  // Proline
constexpr DigitalResidue RES_Q = 13;  // Glutamine
constexpr DigitalResidue RES_R = 14;  // Arginine
constexpr DigitalResidue RES_S = 15;  // Serine
constexpr DigitalResidue RES_T = 16;  // Threonine
constexpr DigitalResidue RES_V = 17;  // Valine
constexpr DigitalResidue RES_W = 18;  // Tryptophan
constexpr DigitalResidue RES_Y = 19;  // Tyrosine

// ============================================================================
// Test Case 1: ConstantAllOnes
// ============================================================================
// Description: All match scores = 1.0, M=5, L=5
// Expected behavior: Each position contributes 1.0, total score depends on
//                    alignment path through DP matrix
// For MSV (ungapped, finding max scoring segment):
// - With all scores = 1.0, the best segment aligns min(M, L) positions
// - Expected score: min(M, L) * 1.0 = 5.0 for full alignment
// ============================================================================

struct ConstantAllOnesTest {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float MATCH_SCORE = 1.0f;
    static constexpr float EXPECTED_SCORE = 5.0f;  // 5 positions * 1.0
    
    // Sequence: A, C, D, E, F (0, 1, 2, 3, 4)
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D, RES_E, RES_F});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 2: ConstantAllTwos
// ============================================================================
// Description: All match scores = 2.0, M=5, L=5
// Expected behavior: Each position contributes 2.0
// Expected score: min(M, L) * 2.0 = 10.0 for full alignment
// ============================================================================

struct ConstantAllTwosTest {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float MATCH_SCORE = 2.0f;
    static constexpr float EXPECTED_SCORE = 10.0f;  // 5 positions * 2.0
    
    // Sequence: G, H, I, K, L (5, 6, 7, 8, 9)
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_G, RES_H, RES_I, RES_K, RES_L});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 3: SinglePositionModel
// ============================================================================
// Description: M=1, L=5 (single position model, longer sequence)
// Model has only one position that matches all residues with score 1.0
// Expected score: 1.0 (only one position to match)
// ============================================================================

struct SinglePositionModelTest {
    static constexpr int MODEL_LENGTH = 1;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float MATCH_SCORE = 1.0f;
    static constexpr float EXPECTED_SCORE = 1.0f;  // Only 1 position in model
    
    // Sequence: A, C, D, E, F
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D, RES_E, RES_F});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 4: SingleResidueSequence
// ============================================================================
// Description: M=5, L=1 (full model, single residue sequence)
// Sequence has only one residue, model has 5 positions
// MSV finds the best matching position for the single residue
// With constant scores = 1.0, expected score: 1.0
// ============================================================================

struct SingleResidueSequenceTest {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 1;
    static constexpr float MATCH_SCORE = 1.0f;
    static constexpr float EXPECTED_SCORE = 1.0f;  // Only 1 residue in sequence
    
    // Sequence: M (Methionine = 10)
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_M});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 5: AlternatingPattern
// ============================================================================
// Description: Position k prefers residue (k % K) with high score
// M=10, L=10
// Sequence is designed to match the pattern: position k gets residue (k-1)%K
// Expected score: sum of match_score for each position = 10 * 3.0 = 30.0
// ============================================================================

struct AlternatingPatternTest {
    static constexpr int MODEL_LENGTH = 10;
    static constexpr int SEQUENCE_LENGTH = 10;
    static constexpr float MATCH_SCORE = 3.0f;
    static constexpr float MISMATCH_SCORE = -1.0f;
    static constexpr float EXPECTED_SCORE = 30.0f;  // 10 positions * 3.0
    
    // Sequence designed to match the pattern: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    // Position k in model prefers residue (k-1), so sequence residue at position i
    // should be (i-1) for best match
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({
            RES_A, RES_C, RES_D, RES_E, RES_F,  // 0, 1, 2, 3, 4
            RES_G, RES_H, RES_I, RES_K, RES_L   // 5, 6, 7, 8, 9
        });
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_alternating_pattern_profile(
            MODEL_LENGTH, MATCH_SCORE, MISMATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 6: AllSameResidue
// ============================================================================
// Description: Sequence "AAAAA" (all Alanine, residue 0)
// M=5, L=5
// Model position 1 gives score 3.0 to A, others give -1.0
// Expected: Best alignment starts at position 1, score = 3.0 + 4*(-1.0) = -1.0
//           Or better: align only position 1 for score 3.0
// For MSV (max segment): best single position = 3.0
// ============================================================================

struct AllSameResidueTest {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 5;
    static constexpr float EXPECTED_SCORE = 3.0f;  // Best single position
    
    // Sequence: A, A, A, A, A (all Alanine = 0)
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_A, RES_A, RES_A, RES_A});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        HMMProfile profile(MODEL_LENGTH, &abc);
        profile.model_length = MODEL_LENGTH;
        profile.name = "all_same_residue_model";
        profile.max_length = 100;
        
        // Only position 1 likes Alanine (score 3.0)
        // All other positions give -1.0 to all residues
        for (int k = 1; k <= MODEL_LENGTH; k++) {
            for (int x = 0; x < abc.K; x++) {
                if (k == 1 && x == RES_A) {
                    profile.match_score(k, x) = 3.0f;  // Position 1 loves Alanine
                } else {
                    profile.match_score(k, x) = -1.0f; // Everything else is bad
                }
            }
        }
        
        return profile;
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 7: AllDifferentResidues
// ============================================================================
// Description: Full alphabet sequence using first 20 residues
// M=20, L=20
// Model and sequence both use full alphabet in order
// Each position k is matched with score 2.0
// Expected score: 20 * 2.0 = 40.0
// ============================================================================

struct AllDifferentResiduesTest {
    static constexpr int MODEL_LENGTH = 20;  // Full alphabet
    static constexpr int SEQUENCE_LENGTH = 20;
    static constexpr float MATCH_SCORE = 2.0f;
    static constexpr float EXPECTED_SCORE = 40.0f;  // 20 positions * 2.0
    
    // Sequence: Full alphabet in order: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({
            RES_A, RES_C, RES_D, RES_E, RES_F, 
            RES_G, RES_H, RES_I, RES_K, RES_L,
            RES_M, RES_N, RES_P, RES_Q, RES_R,
            RES_S, RES_T, RES_V, RES_W, RES_Y
        });
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        HMMProfile profile(MODEL_LENGTH, &abc);
        profile.model_length = MODEL_LENGTH;
        profile.name = "all_different_residues_model";
        profile.max_length = 100;
        
        // Position k prefers residue (k-1) with score MATCH_SCORE
        for (int k = 1; k <= MODEL_LENGTH; k++) {
            int preferred = k - 1;  // Position k prefers residue k-1
            for (int x = 0; x < abc.K; x++) {
                if (x == preferred) {
                    profile.match_score(k, x) = MATCH_SCORE;
                } else {
                    profile.match_score(k, x) = -1.0f;
                }
            }
        }
        
        return profile;
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 8: ShorterSequence
// ============================================================================
// Description: L=3, M=10 (short sequence, longer model)
// Constant scores = 2.0
// MSV can only align 3 positions (sequence limited)
// Expected score: 3 * 2.0 = 6.0
// ============================================================================

struct ShorterSequenceTest {
    static constexpr int MODEL_LENGTH = 10;
    static constexpr int SEQUENCE_LENGTH = 3;
    static constexpr float MATCH_SCORE = 2.0f;
    static constexpr float EXPECTED_SCORE = 6.0f;  // 3 residues * 2.0
    
    // Sequence: A, C, D
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 9: LongerSequence
// ============================================================================
// Description: L=20, M=5 (long sequence, shorter model)
// Constant scores = 1.5
// MSV can only align 5 positions (model limited)
// Expected score: 5 * 1.5 = 7.5
// ============================================================================

struct LongerSequenceTest {
    static constexpr int MODEL_LENGTH = 5;
    static constexpr int SEQUENCE_LENGTH = 20;
    static constexpr float MATCH_SCORE = 1.5f;
    static constexpr float EXPECTED_SCORE = 7.5f;  // 5 positions * 1.5
    
    // Sequence: 20 residues cycling through alphabet
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({
            RES_A, RES_C, RES_D, RES_E, RES_F,
            RES_G, RES_H, RES_I, RES_K, RES_L,
            RES_M, RES_N, RES_P, RES_Q, RES_R,
            RES_S, RES_T, RES_V, RES_W, RES_Y
        });
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        return create_constant_score_profile(MODEL_LENGTH, MATCH_SCORE, abc);
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

// ============================================================================
// Test Case 10: MixedScores (additional test for variety)
// ============================================================================
// Description: M=4, L=4 with specific position-dependent scores
// Position scores:
//   k=1: A=2.0, others=0.0
//   k=2: C=3.0, others=0.0
//   k=3: D=2.0, others=0.0
//   k=4: E=3.0, others=0.0
// Sequence: A, C, D, E
// Expected: 2.0 + 3.0 + 2.0 + 3.0 = 10.0
// ============================================================================

struct MixedScoresTest {
    static constexpr int MODEL_LENGTH = 4;
    static constexpr int SEQUENCE_LENGTH = 4;
    static constexpr float EXPECTED_SCORE = 10.0f;  // 2+3+2+3 = 10.0
    
    // Sequence: A, C, D, E
    static std::vector<DigitalResidue> get_sequence() {
        return create_digital_sequence({RES_A, RES_C, RES_D, RES_E});
    }
    
    static HMMProfile get_profile(const AminoAcidAlphabet& abc) {
        HMMProfile profile(MODEL_LENGTH, &abc);
        profile.model_length = MODEL_LENGTH;
        profile.name = "mixed_scores_model";
        profile.max_length = 100;
        
        // Initialize all scores to 0
        for (int k = 1; k <= MODEL_LENGTH; k++) {
            for (int x = 0; x < abc.K; x++) {
                profile.match_score(k, x) = 0.0f;
            }
        }
        
        // Set specific scores for each position
        profile.match_score(1, RES_A) = 2.0f;  // Position 1 likes A
        profile.match_score(2, RES_C) = 3.0f;  // Position 2 likes C
        profile.match_score(3, RES_D) = 2.0f;  // Position 3 likes D
        profile.match_score(4, RES_E) = 3.0f;  // Position 4 likes E
        
        return profile;
    }
    
    static DPMatrix get_dp_matrix() {
        return DPMatrix(MODEL_LENGTH, SEQUENCE_LENGTH);
    }
};

} // namespace msv_test

#endif // MSV_FILTER_TEST_VECTORS_HPP
