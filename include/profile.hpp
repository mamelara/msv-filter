/*******************************************************************************
 * File: include/profile.hpp
 * Description: HMMER-compatible P7_PROFILE structure for MSV filter mocking.
 * This replicates the essential fields from hmmer/src/hmmer.h
 ******************************************************************************/

#ifndef MSV_FILTER_PROFILE_HPP
#define MSV_FILTER_PROFILE_HPP

#include <vector>
#include <string>
#include <limits>
#include <cstring>
#include <memory>
#include "hmmer_types.hpp"
#include "aa_alphabet.hpp"

/*******************************************************************************
 * P7_PROFILE Structure (from hmmer.h)
 ******************************************************************************/

class HMMProfile {
public:
    // --- Core Data (Pointers replaced with vectors for safety) ---
    
    // Transitions: tsc[0..model_length-1][0..p7P_NTRANS-1]
    // Flattened as: tsc[k * p7P_NTRANS + state]
    std::vector<float> tsc;
    
    // Emissions: rsc[0..Kp-1][0..model_length][0..p7P_NR-1]
    // Access as: rsc[residue_idx][k * p7P_NR + emission_type]
    // MSV only uses MSC (p7P_MSC = 0)
    std::vector<std::vector<float>> rsc;
    
    // Special transitions: xsc[p7P_NXSTATES][p7P_NXTRANS]
    // States: N, E, C, J, B (but E has no transitions in this array)
    // Actually in HMMER: xsc[N,C,J,B][LOOP,MOVE]
    float xsc[p7P_NXSTATES][p7P_NXTRANS];
    
    // --- Configuration & Dimensions ---
    int mode;           // Configured algorithm mode (e.g., p7_LOCAL)
    int sequence_length;              // Current configured target seq length
    int allocM;         // Max # of nodes allocated
    int model_length;              // Actual number of nodes in model
    int max_length;     // Calculated upper bound on emitted seq length
    float nj;           // Expected # of uses of J
    
    // --- Statistical Parameters ---
    float evparam[p7_NEVPARAM];  // E-value parameters
    float cutoff[p7_NCUTOFFS];   // Score cutoffs
    float compo[p7_MAXABET];     // Per-model composition
    
    // --- Alphabet Reference ---
    const AminoAcidAlphabet* abc;
    
    // --- Metadata (simplified for mocking) ---
    std::string name;
    
    // --- Constructor ---
    HMMProfile(int model_length, const AminoAcidAlphabet* alphabet)
        : allocM(model_length), abc(alphabet)
    {
        // Initialize scalars
        mode = 0;
        sequence_length = 0;
        this->model_length = 0;
        max_length = -1;
        nj = 0.0f;
        
        // Initialize special transitions to 0
        for (int i = 0; i < p7P_NXSTATES; i++) {
            for (int j = 0; j < p7P_NXTRANS; j++) {
                xsc[i][j] = 0.0f;
            }
        }
        
        // Initialize statistical parameters
        for (int i = 0; i < p7_NEVPARAM; i++) evparam[i] = 0.0f;
        for (int i = 0; i < p7_NCUTOFFS; i++) cutoff[i] = 0.0f;
        for (int i = 0; i < p7_MAXABET; i++) compo[i] = 0.0f;
        
        // Allocate transitions: allocM * p7P_NTRANS
        // Note: Node 0 has no transitions in HMMER
        tsc.resize(allocM * p7P_NTRANS, -eslINFINITY);
        
        // Allocate emissions: Kp rows, each with (allocM + 1) * p7P_NR columns
        // Note: k=0 has no emissions (nonexistent M_0 and I_0)
        int width = (allocM + 1) * p7P_NR;
        rsc.resize(abc->Kp);
        for (int x = 0; x < abc->Kp; x++) {
            rsc[x].resize(width, -eslINFINITY);
        }
        
        // Set node 0 emissions to -inf (already done by resize)
        // Set gap/missing characters to -inf
        int gap_idx = abc->inmap['-'];
        if (gap_idx < abc->Kp) {
            std::fill(rsc[gap_idx].begin(), rsc[gap_idx].end(), -eslINFINITY);
        }
    }
    
    // --- Accessor Methods (replace HMMER macros) ---
    
    // p7P_TSC(gm, k, s) -> tsc[(k) * p7P_NTRANS + (s)]
    inline float& trans(int k, int state_idx) {
        return tsc[(k * p7P_NTRANS) + state_idx];
    }
    
    inline float trans(int k, int state_idx) const {
        return tsc[(k * p7P_NTRANS) + state_idx];
    }
    
    // p7P_MSC(gm, k, x) -> rsc[x][(k) * p7P_NR + p7P_MSC]
    inline float& match_score(int k, int residue_idx) {
        return rsc[residue_idx][(k * p7P_NR) + p7P_MSC];
    }
    
    inline float match_score(int k, int residue_idx) const {
        return rsc[residue_idx][(k * p7P_NR) + p7P_MSC];
    }
    
    // p7P_ISC(gm, k, x) -> rsc[x][(k) * p7P_NR + p7P_ISC]
    inline float& insert_score(int k, int residue_idx) {
        return rsc[residue_idx][(k * p7P_NR) + p7P_ISC];
    }
};

#endif // MSV_FILTER_PROFILE_HPP
