/*******************************************************************************
 * File: include/mock_data.hpp
 * Description: Mock data generator for MSV filter testing.
 * Creates hardcoded test inputs for p7_GMSV function.
 ******************************************************************************/

#ifndef MSV_FILTER_MOCK_DATA_HPP
#define MSV_FILTER_MOCK_DATA_HPP

#include <vector>
#include <string>
#include <cmath>
#include "hmmer_types.hpp"
#include "aa_alphabet.hpp"
#include "profile.hpp"
#include "dp_matrix.hpp"

/*******************************************************************************
 * Mock Data Generator
 * 
 * Provides hardcoded test data for:
 * - Digital sequences (DigitalResidue arrays)
 * - Profiles with match scores
 * - DP matrices
 ******************************************************************************/

class MockDataGenerator {
public:
    // --- Mock Digital Sequence ---
    // Creates a simple 1-indexed digital sequence
    // Returns vector where index 0 is sentinel, 1..sequence_length are residues
    static std::vector<DigitalResidue> create_simple_sequence(int sequence_length, const AminoAcidAlphabet& abc) {
        std::vector<DigitalResidue> digital_sequence(sequence_length + 2);  // 0 = sentinel, 1..sequence_length = seq, sequence_length+1 = sentinel
        digital_sequence[0] = digitalResidueSentinel;
        digital_sequence[sequence_length + 1] = digitalResidueSentinel;
        
        // Simple pattern: A(0), C(1), D(2), E(3), ... wrapping around
        for (int i = 1; i <= sequence_length; i++) {
            digital_sequence[i] = (i - 1) % abc.K;
        }
        
        return digital_sequence;
    }
    
    // --- Mock Profile ---
    // Creates a profile with model_length nodes and simple match scores
    // MSV only uses match scores (MSC), not transitions
    static HMMProfile create_simple_profile(int model_length, const AminoAcidAlphabet& abc) {
        HMMProfile profile(model_length, &abc);
        profile.model_length = model_length;
        profile.name = "test_model";
        
        // Set max_length (used in transitions calculations)
        profile.max_length = 100;
        
        // Fill match scores with simple pattern
        // For each position k (1..model_length) and each residue x (0..K-1):
        // score = sin(k + x) * 2.0  (range: -2 to 2)
        for (int k = 1; k <= model_length; k++) {
            for (int x = 0; x < abc.K; x++) {
                float score = std::sin(float(k + x)) * 2.0f;
                profile.match_score(k, x) = score;
            }
        }
        
        return profile;
    }
    
    // --- Mock Profile with Constant Scores ---
    // Creates a profile where all match scores are the same value
    static HMMProfile create_constant_profile(int model_length, const AminoAcidAlphabet& abc, float score) {
        HMMProfile profile(model_length, &abc);
        profile.model_length = model_length;
        profile.name = "constant_model";
        profile.max_length = 100;
        
        for (int k = 1; k <= model_length; k++) {
            for (int x = 0; x < abc.K; x++) {
                profile.match_score(k, x) = score;
            }
        }
        
        return profile;
    }
    
    // --- Mock Profile with Specific Pattern ---
    // Creates a profile with a specific residue pattern at each position
    // Useful for testing specific alignment scenarios
    static HMMProfile create_pattern_profile(int model_length, const AminoAcidAlphabet& abc) {
        HMMProfile profile(model_length, &abc);
        profile.model_length = model_length;
        profile.name = "pattern_model";
        profile.max_length = 100;
        
        // Pattern: position k prefers residue (k % K)
        for (int k = 1; k <= model_length; k++) {
            int preferred = (k - 1) % abc.K;
            for (int x = 0; x < abc.K; x++) {
                if (x == preferred) {
                    profile.match_score(k, x) = 2.0f;  // High score for match
                } else {
                    profile.match_score(k, x) = -1.0f;  // Penalty for mismatch
                }
            }
        }
        
        return profile;
    }
    
    // --- Create DP Matrix ---
    static DPMatrix create_dp_matrix(int model_length, int sequence_length) {
        return DPMatrix(model_length, sequence_length);
    }
    
    // --- Create Test Case: Simple ---
    // sequence_length=10, model_length=5, simple sequence and profile
    static void create_simple_test(
        std::vector<DigitalResidue>& digital_sequence,
        int& sequence_length,
        HMMProfile& profile,
        DPMatrix& dp_matrix,
        const AminoAcidAlphabet& abc
    ) {
        sequence_length = 10;
        int model_length = 5;
        
        digital_sequence = create_simple_sequence(sequence_length, abc);
        profile = create_simple_profile(model_length, abc);
        dp_matrix = create_dp_matrix(model_length, sequence_length);
    }
    
    // --- Create Test Case: Constant ---
    // All match scores = 1.0, should give predictable results
    static void create_constant_test(
        std::vector<DigitalResidue>& digital_sequence,
        int& sequence_length,
        HMMProfile& profile,
        DPMatrix& dp_matrix,
        const AminoAcidAlphabet& abc
    ) {
        sequence_length = 20;
        int model_length = 10;
        
        digital_sequence = create_simple_sequence(sequence_length, abc);
        profile = create_constant_profile(model_length, abc, 1.0f);
        dp_matrix = create_dp_matrix(model_length, sequence_length);
    }
    
    // --- Create Test Case: Pattern ---
    // Sequence and profile with matching patterns
    static void create_pattern_test(
        std::vector<DigitalResidue>& digital_sequence,
        int& sequence_length,
        HMMProfile& profile,
        DPMatrix& dp_matrix,
        const AminoAcidAlphabet& abc
    ) {
        sequence_length = 15;
        int model_length = 10;
        
        digital_sequence = create_simple_sequence(sequence_length, abc);
        profile = create_pattern_profile(model_length, abc);
        dp_matrix = create_dp_matrix(model_length, sequence_length);
    }
    
    // --- Utility: Print Digital Sequence ---
    static void print_sequence(const std::vector<DigitalResidue>& digital_sequence, int sequence_length, const AminoAcidAlphabet& abc) {
        std::cout << "Digital Sequence (length " << sequence_length << "): ";
        for (int i = 1; i <= sequence_length; i++) {
            int residue = digital_sequence[i];
            if (residue < abc.K) {
                std::cout << abc.sym[residue];
            } else {
                std::cout << "?";
            }
        }
        std::cout << std::endl;
    }
    
    // --- Utility: Print Profile Match Scores ---
    static void print_profile(const HMMProfile& profile, int max_k = -1) {
        int model_length = (max_k > 0) ? std::min(max_k, profile.model_length) : profile.model_length;
        std::cout << "Profile (model_length=" << profile.model_length << "):" << std::endl;
        std::cout << "  First 5 match scores at each position:" << std::endl;
        for (int k = 1; k <= model_length; k++) {
            std::cout << "    k=" << k << ": ";
            for (int x = 0; x < std::min(5, 20); x++) {
                std::cout << profile.match_score(k, x) << " ";
            }
            std::cout << "..." << std::endl;
        }
    }
};

#endif // MSV_FILTER_MOCK_DATA_HPP
