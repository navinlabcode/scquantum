#include <stdint.h>
#include <stdio.h>
#include <tgmath.h>
#define FELD_PI 3.141592653589793
int ecf_local_max(double *s, int *nys, double *ys, int *nsigmas, double *sigmas, double *peak_location, double *peak_height, double *peak_phase)
{
    double v0;
    uint32_t v1;
    uint32_t v3;
    uint32_t v5;
    uint32_t v7;
    uint32_t let9;
    double let10;
    double let11;
    double let12;
    double state13;
    uint32_t v14;
    double r91;
    double let92;
    double let93;
    double _Complex state94;
    uint32_t v95;
    double state97;
    uint32_t v98;
    double _Complex r100;
    
    v0 = *s;
    /* printf("s = %.10f\n", v0); */

    v1 = *nys;
    /* printf("Number of observations = %d\n", v1); */
    
    double _a2[v1];
    double *a2 = _a2;
    
    a2 = ys;
    /* for (v3 = 0; v3 < v1; v3++) {
        printf("Observation number %d: %.10f\n", v3, a2[v3]);
    } */
    v5 = *nsigmas;
    /* printf("Number of standard deviations: %d\n", v5); */
    
    double _a6[v5];
    double *a6 = _a6;
    
    a6 = sigmas;
    /* for (v7 = 0; v7 < v5; v7++) {
        printf("Standard deviation number %d: %.10f\n", v7, a6[v7]);
    } */
    let9 = v5 <= v1 ? v5 : v1;
    let10 = 2.0 * FELD_PI;
    let11 = 4.0 * FELD_PI;
    let12 = FELD_PI * FELD_PI;
    state13 = v0;
    for (v14 = 0; v14 < 20; v14++) {
        double let15;
        double let16;
        double _Complex state17;
        uint32_t v18;
        double _Complex let24;
        double let25;
        double let26;
        double let27;
        double _Complex state28;
        uint32_t v29;
        double _Complex let41;
        double state42;
        uint32_t v43;
        double let45;
        double state46;
        uint32_t v47;
        double let52;
        double let53;
        double state54;
        uint32_t v55;
        double let62;
        double let63;
        double _Complex state64;
        uint32_t v65;
        double _Complex let67;
        double _Complex let68;
        double _Complex let69;
        double _Complex let70;
        double _Complex let71;
        double let72;
        double _Complex let73;
        double _Complex let74;
        double _Complex let75;
        double _Complex let76;
        double _Complex let77;
        double _Complex let78;
        double let79;
        double _Complex let80;
        double _Complex let81;
        double _Complex let82;
        double let83;
        double _Complex let84;
        double _Complex let85;
        double _Complex let86;
        double _Complex let87;
        double let88;
        double _Complex let89;
        double _Complex let90;
        
        let15 = let10 * state13;
        let16 = let11 * state13;
        state17 = 0.0;
        for (v18 = 0; v18 < let9; v18++) {
            double let19;
            double let20;
            double let21;
            double let22;
            double let23;
            
            let19 = let15 * a6[v18];
            let20 = exp(let19 * let19);
            let21 = a2[v18];
            let22 = -1.0 + let20;
            let23 = a6[v18];
            state17 = 6.283185307179586 * I * (let20 * exp(I * (let15 *
                                                                let21))) *
                (let22 * let21 + I * (let16 * (let23 * let23))) * (1.0 /
                                                                   (let22 *
                                                                    let22)) +
                state17;
        }
        let24 = state17;
        let25 = FELD_PI * state13;
        let26 = let11 * state13;
        let27 = 78.95683520871486 * (state13 * state13);
        state28 = 0.0;
        for (v29 = 0; v29 < let9; v29++) {
            double let30;
            double let31;
            double let32;
            double let33;
            double let34;
            double let35;
            double let36;
            double let37;
            double let38;
            double _Complex let39;
            double let40;
            
            let30 = let15 * a6[v29];
            let31 = -1.0 + exp(let30 * let30);
            let32 = let15 * a6[v29];
            let33 = exp(let32 * let32);
            let34 = a2[v29];
            let35 = let25 * a6[v29];
            let36 = let34 * let34;
            let37 = a6[v29];
            let38 = let37 * let37;
            let39 = let34 + I * -(let26 * let38);
            let40 = a6[v29];
            state28 = -(39.47841760435743 / (let31 * let31 * let31)) * (let33 *
                                                                        exp(I *
                                                                        (let15 *
                                                                         let34))) *
                (exp(8.0 * (let35 * let35)) * let36 - 2.0 * let38 + let39 *
                 let39 - 2.0 * let33 * (let36 - let38 + let27 * (let38 * let40 *
                                                                 let40) + I *
                                        -(let26 * let34 * let38))) + state28;
        }
        let41 = state28;
        state42 = 0.0;
        for (v43 = 0; v43 < v5; v43++) {
            double let44;
            
            let44 = let10 * a6[v43] * state13;
            state42 = 1.0 / (1.0 - exp(-(let44 * let44))) + state42;
        }
        let45 = state42;
        state46 = 0.0;
        for (v47 = 0; v47 < v5; v47++) {
            double let48;
            double let49;
            double let50;
            double let51;
            
            let48 = let15 * a6[v47];
            let49 = exp(-(let48 * let48));
            let50 = a6[v47];
            let51 = 1.0 - let49;
            state46 = -(8.0 * let49 * let12 * state13 * (let50 * let50) /
                        (let51 * let51)) + state46;
        }
        let52 = state46;
        let53 = FELD_PI * state13;
        state54 = 0.0;
        for (v55 = 0; v55 < v5; v55++) {
            double let56;
            double let57;
            double let58;
            double let59;
            double let60;
            double let61;
            
            let56 = let15 * a6[v55];
            let57 = exp(let56 * let56);
            let58 = a6[v55] * FELD_PI;
            let59 = let53 * a6[v55];
            let60 = 8.0 * (let59 * let59);
            let61 = -1.0 + let57;
            state54 = 8.0 * let57 * (let58 * let58) * (1.0 + let60 + let57 *
                                                       (-1.0 + let60)) /
                (let61 * let61 * let61) + state54;
        }
        let62 = state54;
        let63 = 0.5 * state13;
        state64 = 0.0;
        for (v65 = 0; v65 < let9; v65++) {
            double let66;
            
            let66 = let10 * a6[v65] * state13;
            state64 = 1.0 / (1.0 - exp(-(let66 * let66))) * exp(I * (let15 *
                                                                     a2[v65])) +
                state64;
        }
        let67 = state64;
        let68 = 2.0 * let67;
        let69 = let24;
        let70 = 2.0 * let69;
        let71 = let41;
        let72 = let45;
        let73 = let67 / let72;
        let74 = let72 * let69;
        let75 = 1.0 / (let72 * let72);
        let76 = let72 * let72 * let71;
        let77 = let72;
        let78 = 1.0 / (let72 * let72 * let72);
        let79 = let52;
        let80 = (let74 - let67 * let79) * let75;
        let81 = let68 * (let79 * let79) + let76;
        let82 = let70 * let79;
        let83 = let62;
        let84 = (let81 - let77 * (let82 + let67 * let83)) * let78;
        let85 = let73;
        let86 = conj(let85);
        let87 = let80;
        let88 = creal(conj(let87) * let85 + let86 * let87);
        let89 = 2.0 * conj(let87) * let87;
        let90 = let84;
        state13 = let63 + 0.5 * (state13 - let88 / creal(conj(let90) * let85 +
                                 let89 + let86 * let90));
    }
    r91 = state13;
    let92 = 2.0 * FELD_PI;
    let93 = let92 * r91;
    state94 = 0.0;
    for (v95 = 0; v95 < (v5 <= v1 ? v5 : v1); v95++) {
        double let96;
        
        let96 = let92 * a6[v95] * r91;
        state94 = 1.0 / (1.0 - exp(-(let96 * let96))) * exp(I * (let93 *
                                                                 a2[v95])) +
            state94;
    }
    state97 = 0.0;
    for (v98 = 0; v98 < v5; v98++) {
        double let99;
        
        let99 = let92 * a6[v98] * r91;
        state97 = 1.0 / (1.0 - exp(-(let99 * let99))) + state97;
    }
    r100 = state94 / state97;
    *peak_location = r91; 
    *peak_height = cabs(r100);
    *peak_phase = carg(r100);
    /* printf("Peak location is %.3f, height is %.3f, phase is %.3f\n", *peak_location, *peak_height, *peak_phase); */
    return 0;
}

