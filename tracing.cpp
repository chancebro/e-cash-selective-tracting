// tracing.cpp

#include "tracing.h"
#include <functional>

using namespace mcl::bn;
using namespace std::chrono;

mcl::bn::Fr alpha;
mcl::bn::Fr beta;
mcl::bn::G2 W;
mcl::bn::G2 X;
sqlite3* db = nullptr;  // 전역 변수 정의

void hashToG1(mcl::bn::G1& g, const std::string& label) {
    uint32_t seed = std::hash<std::string>{}(label);
    mapToG1(g, seed);
}

void hashToG2(mcl::bn::G2& h, const std::string& label) {
    uint32_t seed = std::hash<std::string>{}(label);
    mapToG2(h, seed);
}

void hashToFr(mcl::bn::Fr &target,const std::string& msg){
    target.setHashOf(msg);  // MCL 내장 함수
}
void hashToGT(mcl::bn::Fp12& target,const std::string& msg) {

    G1 G1;
    hashToG1(G1,msg);
    G2 G2; 
    hashToG2(G2,msg + "#");  // 도메인 분리
    pairing(target, G1, G2);
    
}
void setupG_Generators(G1& g, G1& g0, G1& g1, G1& ge,G2& gk) {
    hashToG1(g,  "g");
    hashToG1(g0, "g0");
    hashToG1(g1, "g1");
    hashToG1(ge, "ge");
    hashToG2(gk, "gk");

}

void setupH_Generators(G1& h, G1& h0, G1& h1, G1& h2,G1 &he,G1 &ht,G2 &g) {
    hashToG1(h,  "h");
    hashToG1(h0, "h0");
    hashToG1(h1, "h1");
    hashToG1(h2, "h2");
    hashToG1(he, "he");
    hashToG1(ht, "ht");
    hashToG2(g,"g");
}

// Prover: 증명 생성
pk_1_proof pk_1_prove(const Fr &sp, const Fr &u, const G1 &g0, const G1 &g1, const G1 &Cm,const Fp12 &G) {
    pk_1_proof proof;
    Fr r_sp, r_u;
    r_sp.setByCSPRNG();
    r_u.setByCSPRNG();

    G1::mul(proof.Cm_sp, g0, r_sp);
    G1::mul(proof.Cm_u, g1, r_u);
    Fp12::pow(proof.U, G,u);
    Fp12::pow(proof.Cm_U, G,r_u);

    proof.c.setByCSPRNG();
    proof.Z_sp = r_sp + proof.c * sp;
    proof.Z_u  = r_u  + proof.c * u;

    return proof;
}

// verify pk1 for account est
bool pk_1_verify(const pk_1_proof &proof, const G1 &Cm, const G1 &g0, const G1 &g1,const Fp12 &G) {
    G1 lhs_g0, lhs_g1, lhs;
    G1::mul(lhs_g0, g0, proof.Z_sp);
    G1::mul(lhs_g1, g1, proof.Z_u);
    G1::add(lhs, lhs_g0, lhs_g1);

    G1 rhs_Cm, rhs_tmp, rhs;
    G1::mul(rhs_Cm, Cm, proof.c);
    G1::add(rhs_tmp, proof.Cm_sp, proof.Cm_u);
    G1::add(rhs, rhs_tmp, rhs_Cm);

    Fp12 U_c,lU;
    Fp12::pow(U_c,proof.U,proof.c);
    Fp12::mul(lU, U_c, proof.Cm_U);
    Fp12 rU;
    Fp12::pow(rU,G,proof.Z_u);

    return (lhs == rhs)&&(lU==rU);
}
bool account_est(User &user_out, mcl::bn::G1& g, mcl::bn::G1& g0, mcl::bn::G1& g1, mcl::bn::G2& g_G2,const mcl::bn::Fp12& G) {
    // Step 1: User generates commitment Cm
    //user side 1 start
    Fr sp;
    sp.setByCSPRNG();
    

    G1 g0sp, g1u, Cm;
    G1::mul(g0sp, g0, sp);
    G1::mul(g1u, g1, user_out.u);
    G1::add(Cm, g0sp, g1u);

    // Step 2: User proves correctness of Cm
    pk_1_proof proof = pk_1_prove(sp, user_out.u, g0, g1, Cm,G);
    if (!pk_1_verify(proof, Cm,g0, g1,G)) {
        std::cout << "❌ Verification failed: Pk_1 is not valid" << std::endl;
        return false;
    }
    //user side 1 end

    // Step 3: Bank sign on Cm for gen account
    //bank side start
    Fr spp, e;
    spp.setByCSPRNG();
    e.setByCSPRNG();

    G1 g0spp, gg0spp, Cm_gg0spp;
    G1::mul(g0spp, g0, spp);
    G1::add(gg0spp, g, g0spp);
    G1::add(Cm_gg0spp, Cm, gg0spp);

    Fr denom = e + alpha;
    Fr inv_denom;
    Fr::inv(inv_denom, denom);

    G1 A;
    G1::mul(A, Cm_gg0spp, inv_denom);//A is signature of bank  
    //bank side end

    //user side 2 start
    // Step 4: Store into User struct
    user_out.A = A;
    user_out.e = e;
    user_out.s = sp + spp;

    // Step 5: Signature verification
    G2 W, he, Whe;
    G2::mul(W, g_G2, alpha);
    G2::mul(he, g_G2, e);
    G2::add(Whe, he, W);

    G1 g0s, gg0s, gg0sg1u;
    G1::mul(g0s, g0, user_out.s);
    G1::add(gg0s, g, g0s);
    G1::add(gg0sg1u, gg0s, g1u);

    Fp12 e1, e2;
    pairing(e1, user_out.A, Whe);
    pairing(e2, gg0sg1u, g_G2);

    if (e1 == e2) {
        //std::cout << "✅ Account establishment successful!" << std::endl;
        return true;
    } else {
        //std::cout << "❌ Signature pairing verification failed!" << std::endl;
        return false;
    }
    //user side 2 end
}
// Prover: 증명 생성
pk_2_proof pk_2_prove(const mcl::bn::Fr &tp, const mcl::bn::Fr &u,const mcl::bn::Fr &vp,
    const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2, const mcl::bn::G1 &Cmp,const mcl::bn::Fp12 &G) {
    pk_2_proof proof;

    Fr r_tp, r_u,r_vp;
    r_tp.setByCSPRNG();
    r_u.setByCSPRNG();
    r_vp.setByCSPRNG();

    G1::mul(proof.Cm_tp, h0, r_tp);
    G1::mul(proof.Cm_u, h1, r_u);
    G1::mul(proof.Cm_vp,h2,r_vp);
    Fp12::pow(proof.Cm_U,G,r_u);
    Fp12::pow(proof.U,G,u);

    proof.c.setByCSPRNG();
    proof.Z_tp = r_tp + proof.c * tp;
    proof.Z_u  = r_u  + proof.c * u;
    proof.Z_vp = r_vp + proof.c * vp;
    proof.Z_U = r_u+proof.c*u;

    return proof;
}

// Verifier: 증명 검증
bool pk_2_verify(const pk_2_proof &proof, const mcl::bn::G1 &Cmp, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1 , const mcl::bn::G1 &h2 ,const mcl::bn::Fp12 &G) {

    G1 lhs_h0, lhs_h1, lhs_h2,temp, lhs;
    G1::mul(lhs_h0, h0, proof.Z_tp);
    G1::mul(lhs_h1, h1, proof.Z_u);
    G1::mul(lhs_h2, h2, proof.Z_vp);
    G1::add(temp,lhs_h0,lhs_h1);
    G1::add(lhs,temp,lhs_h2);

    G1 rhs_Cmp, rhs_temp1, rhs_temp2,rhs;
    G1::mul(rhs_Cmp, Cmp, proof.c);
    G1::add(rhs_temp1, proof.Cm_tp, proof.Cm_u);
    G1::add(rhs_temp2, rhs_temp1, rhs_Cmp);
    G1::add(rhs,rhs_temp2,proof.Cm_vp);

    Fp12 U_c,lU;
    Fp12::pow(U_c,proof.U,proof.c);
    Fp12::mul(lU, U_c, proof.Cm_U);

    Fp12 rU1,rU2,rU;
    Fp12::pow(rU,G,proof.Z_u);
    
    return (lhs == rhs)&&(lU==rU) ;
}
bool withdraw_coin(User &user_out, mcl::bn::G1& h, mcl::bn::G1& h0, mcl::bn::G1& h1,mcl::bn::G1& h2, mcl::bn::G2& h_G2,const mcl::bn::Fp12& G) {
    //user side 1 start
    // Step 1: User generates commitment Cm
    //user pick t' v' randomly
    Fr tp,vp;
    tp.setByCSPRNG();
    vp.setByCSPRNG();
    G1 h0_tp,h1_u,h2_vp,temp,Cmp;
    //make Cm' 
    G1::mul(h0_tp,h0,tp);
    G1::mul(h1_u,h1,user_out.u);
    G1::mul(h2_vp,h2,vp);
    G1::add(temp,h0_tp,h1_u);
    G1::add(Cmp,h2_vp,temp);

    //proof Cmp and U
    pk_2_proof proof = pk_2_prove(tp, user_out.u, vp, h0, h1,h2, Cmp,G);
    if (!pk_2_verify(proof, Cmp, h0, h1,h2,G)) {
        std::cout << "❌ Verification failed: Pk_2 is not valid" << std::endl;
        return false;
    }
    //user side 1 end
    //bank side start
    //bank pick for sign on coin
    G1 h0_tpp,h2_vpp;
    Fr tpp,vpp,f;
    tpp.setByCSPRNG();
    vpp.setByCSPRNG();
    f.setByCSPRNG();
    G1::mul(h0_tpp,h0,tpp);//h0^t''
    G1::mul(h2_vpp,h2,vpp);//h2^v''
    G1::add(temp,h,h0_tpp);//h*h0^t''
    G1::add(temp,h2_vpp,temp);//h*h0^t''*h2^v''
    G1::add(temp,Cmp,temp);//h*Cm'*h0^t''*h2^v''
    Fr denom = f + beta;//beta+f
    Fr inv_denom;
    Fr::inv(inv_denom, denom);//(1/beta+1)
    G1::mul(user_out.B,temp,inv_denom);//sign

    //bank side end

    //usuer side 2 start
    //user verify generated signature
    Fr t,v;
    user_out.f=f;
    user_out.t=tp+tpp;
    user_out.v=vp+vpp;
    G2 X, hf, Xhf;//X is bank's pubkey
    G2::mul(X, h_G2, beta);
    G2::mul(hf, h_G2, f);
    G2::add(Xhf, hf, X);
    
    //user verify vaidation of coin
    G1 h0t,h2v, h0th1u,h0th1uh2v,hh0th1uh2v;
    G1::mul(h0t, h0, user_out.t);
    G1::mul(h2v, h2, user_out.v);
    G1::add(h0th1u, h0t, h1_u);
    G1::add(h0th1uh2v, h0th1u, h2v);
    G1::add(hh0th1uh2v,h0th1uh2v, h);
    //verifing 
    Fp12 e1, e2;
    pairing(e1, user_out.B, Xhf);
    pairing(e2, hh0th1uh2v, h_G2);

    if (e1 == e2) {
        //std::cout << "✅ Coin withdraw successful!" << std::endl;
        return true;
    } else {
        //std::cout << "❌ Coin withdraw  failed!" << std::endl;
        return false;
    }
    //user side 2 end
}
pk_3_proof pk_3_prove(User &user_out,const mcl::bn::G1 g, const mcl::bn::G1 &g0, const mcl::bn::G1 g1,const mcl::bn::G2 g_G2, const mcl::bn::Fp12 &G,std::string INFO) {
    Fr r1,r2;
    Fr delta1,delta2;
    G1 g0r1,g1r2,g1r1;
    G1 g0d1,g1d2;//d1 mean delta1,d2 mean delta2
    
    G2 W;
    G2::mul(W,g_G2,alpha);

    r1.setByCSPRNG();
    r2.setByCSPRNG();
    delta1=r1*user_out.e;
    delta2=r2*user_out.e;

    pk_3_proof proof;
    //user n하고 info저장
    user_out.N.setByCSPRNG();
    proof.N=user_out.N;
    //proof를 위한 n하고 ifo저장
    user_out.INFO=INFO;
    proof.INFO=INFO;

    std::string input = user_out.INFO + user_out.N.getStr(10);
    Fp12 G_t_info;
    hashToGT(G_t_info,input);

    //M생성 후 proof와 user에 저장장
    Fp12::pow(proof.M,G_t_info,user_out.u);
    Fp12::pow(user_out.M,G_t_info,user_out.u);

    G1::mul(g0r1,g0,r1);
    G1::mul(g1r2,g1,r2);
    G1::mul(g1r1,g1,r1);

    G1::add(proof.A1,g0r1,g1r2);
    G1::mul(proof.A1e,proof.A1,user_out.e);
    G1::add(proof.A2,user_out.A,g1r1);
    
    //pick random value for pk3
    Fr r_r1,r_r2,r_d1,r_d2,r_s,r_u,r_e,c;

    r_r1.setByCSPRNG();
    r_r2.setByCSPRNG();
    r_d1.setByCSPRNG();
    r_d2.setByCSPRNG();
    r_s.setByCSPRNG();
    r_u.setByCSPRNG();
    r_e.setByCSPRNG();
    c.setByCSPRNG();
    proof.c=c;
   
    //gen Commit
    //G1 g0r1,g1r2;
    Fp12 ee,er1,ed1,es,eu;

    G1 g0r_r1,g1r_r2;
    G1::mul(g0r_r1,g0,r_r1);
    G1::mul(g1r_r2,g1,r_r2);
    G1::add(proof.Cm_A1,g0r_r1,g1r_r2);

    G1::mul(proof.Cm_A1e,proof.A1,r_e);

    pairing(ee,proof.A2,g_G2);
    Fp12::pow(proof.Cm_ee,ee,r_e);
    pairing(er1,g1,W);
    Fp12::pow(proof.Cm_er1,er1,r_r1);
    pairing(ed1,g1,g_G2);
    Fp12::pow(proof.Cm_ed1,ed1,r_d1);
    pairing(es,g0,g_G2);
    Fp12::pow(proof.Cm_es,es,r_s);
    pairing(eu,g1,g_G2);
    Fp12::pow(proof.Cm_eu,eu,r_u);


    Fp12::pow(proof.Cm_M,G_t_info,r_u);

    //gen Z
    proof.Z_r1=r_r1+c*r1;
    proof.Z_r2=r_r2+c*r2;
    proof.Z_d1=r_d1+c*delta1;
    proof.Z_d2=r_d2+c*delta2;
    proof.Z_s=r_s+c*user_out.s;
    proof.Z_u=r_u+c*user_out.u;
    proof.Z_e=r_e+c*user_out.e;

    return proof;
}
bool pk_3_verify(const pk_3_proof &proof, const mcl::bn::G1 g, const mcl::bn::G1 &g0, const mcl::bn::G1 g1,const mcl::bn::G2 g_G2, const mcl::bn::Fp12 &G) {

    G1 lA1,rA1,A1c,g0Z_r1,g1Z_r2;
    G1::mul(A1c,proof.A1,proof.c);
    G1::add(lA1,proof.Cm_A1,A1c);
    G1::mul(g0Z_r1,g0,proof.Z_r1);
    G1::mul(g1Z_r2,g1,proof.Z_r2);
    G1::add(rA1,g0Z_r1,g1Z_r2);
    //std::cout << "[A1 검증] " << (lA1 == rA1 ? "✅ 성공" : "❌ 실패") << std::endl;

    G1 lA1e,rA1e,A1ec,g0Z_d1,g1Z_d2;
    G1::mul(A1ec,proof.A1e,proof.c);
    G1::add(lA1e,proof.Cm_A1e,A1ec);
    G1::mul(rA1e,proof.A1,proof.Z_e);

    //G1::mul(g1Z_d2,g1,proof.Z_d2); 
    //G1::add(rA1e,g0Z_d1,g1Z_d2);

    //std::cout << "[A1e 검증] " << (lA1e == rA1e ? "✅ 성공" : "❌ 실패") << std::endl;

    Fp12 e_A2w, e_gg_G2, result;

    pairing(e_A2w, proof.A2, W);       // e(A2, w)
    pairing(e_gg_G2, g, g_G2);     // e(g0, h0)

    Fp12 inv_e_gg_G2;
    Fp12::inv(inv_e_gg_G2, e_gg_G2);   // 역원 계산
    Fp12::mul(result, e_A2w, inv_e_gg_G2); // 최종 결과 계산

    Fp12 templ1,templ2,templ3,templ4,el,resultc,inv_Cm_ee;

    Fp12::pow(resultc,result,proof.c);

    Fp12::inv(inv_Cm_ee,proof.Cm_ee);

    Fp12::mul(templ1,inv_Cm_ee,proof.Cm_er1);

    Fp12::mul(templ2,proof.Cm_ed1,templ1);
    Fp12::mul(templ3,proof.Cm_es,templ2);
    Fp12::mul(templ4,proof.Cm_eu,templ3);
    Fp12::mul(el,templ4,resultc);

    Fp12 tempr1,tempr2,tempr3,tempr4,tempr5,inv_tempr1,er;

    Fp12 ee,er1,ed1,es,eu;

    pairing(ee,proof.A2,g_G2);
    pairing(er1,g1,W);
    pairing(ed1,g1,g_G2);
    pairing(es,g0,g_G2);
    pairing(eu,g1,g_G2);


    Fp12::pow(tempr1,ee,proof.Z_e);
    Fp12::inv(inv_tempr1, tempr1);   // 역원 계산
    Fp12::pow(tempr2,er1,proof.Z_r1);
    Fp12::pow(tempr3,ed1,proof.Z_d1);
    Fp12::pow(tempr4,es,proof.Z_s);
    Fp12::pow(tempr5,eu,proof.Z_u);
    Fp12 t1,t2,t3;
    Fp12::mul(t1,inv_tempr1,tempr2);
    Fp12::mul(t2,t1,tempr3);
    Fp12::mul(t3,t2,tempr4);
    Fp12::mul(er,t3,tempr5);

    //std::cout << "[페어링 검증] " << (er == el ? "✅ 성공" : "❌ 실패") << std::endl;

    //gen G(Info||N)


    std::string input = proof.INFO + proof.N.getStr(10);
    Fp12 G_t_info;
    hashToGT(G_t_info,input);

    Fp12 M_c,lM;
    Fp12::pow(M_c,proof.M,proof.c);
    Fp12::mul(lM, M_c, proof.Cm_M);

    Fp12 rM;
    Fp12::pow(rM,G_t_info ,proof.Z_u);

    //std::cout << "[M 검증] " << (lU==rU ? "✅ 성공" : "❌ 실패") << std::endl;


    return (lA1 == rA1)&&(lA1e == rA1e)&&(er == el)&&(lM==rM) ;
}


SPK_proof SPK_prove(User &user_out,const mcl::bn::G1 h, const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2,const mcl::bn::G2 h_G2, const mcl::bn::Fp12 &G,const mcl::bn::Fp12 &H,const mcl::bn::Fp12 &H1,std::string INFO,pk_3_proof pk_3) {
        
    Fr r1,r2;
    Fr delta1,delta2;
    G1 h0r1,h1r2,h1r1;
    G1 h0d1,h1d2;//d1 mean delta1,d2 mean delta2
    
    //pick random r1,r1 and calculate delta1,delta2 for spk

    r1.setByCSPRNG();
    r2.setByCSPRNG();
    delta1=r1*user_out.f;
    delta2=r2*user_out.f;

    SPK_proof proof;
    

    //calculate S,D

    //gen S
    Fp12::pow(proof.S,H,user_out.v);
    //gen D
    Fp12 temp_D1,temp_D2;
    std::string input = INFO + pk_3.N.getStr(10)+pk_3.M.getStr(10);
    proof.N=pk_3.N;
    proof.M=pk_3.M;
    Fr R;
    hashToFr(R,input);
    proof.R=R;
    Fp12::pow(temp_D1,H1,R);
    Fp12::pow(temp_D2,temp_D1,user_out.v);
    Fp12::mul(proof.D,user_out.U,temp_D2);

    
    

    

    G1::mul(h0r1,h0,r1);
    G1::mul(h1r2,h1,r2);
    G1::mul(h1r1,h1,r1);

    G1::add(proof.B1,h0r1,h1r2);
    G1::mul(proof.B1f,proof.B1,user_out.f);
    G1::add(proof.B2,user_out.B,h1r1);
    
    //pick random value for spk
    Fr r_r1,r_r2,r_d1,r_d2,r_t,r_u,r_v,r_f,c;

    r_r1.setByCSPRNG();
    r_r2.setByCSPRNG();
    r_d1.setByCSPRNG();
    r_d2.setByCSPRNG();
    r_t.setByCSPRNG();
    r_u.setByCSPRNG();
    r_f.setByCSPRNG();
    r_v.setByCSPRNG();
    c.setByCSPRNG();
    proof.c=c;
    //gen Cm_S and Cm_D
    Fp12 temp_cmD1,temp_cmD2;
    //gen Cm_S
    Fp12::pow(proof.Cm_S,H,r_v);
    //gen Cm_D
    Fp12::pow(temp_cmD1,G,r_u);
    Fp12::pow(temp_cmD2,H1,R);
    Fp12::pow(temp_cmD2,temp_cmD2,r_v);
    Fp12::mul(proof.Cm_D,temp_cmD1,temp_cmD2);
   
    //gen Commit for prove pairing
    //G1 g0r1,g1r2;
    Fp12 ef,er1,ed1,et,eu,ev;

    G1 h0r_r1,h1r_r2;
    G1::mul(h0r_r1,h0,r_r1);
    G1::mul(h1r_r2,h1,r_r2);
    G1::add(proof.Cm_B1,h0r_r1,h1r_r2);

    G1::mul(proof.Cm_B1f,proof.B1,r_f);


    pairing(ef,proof.B2,h_G2);
    Fp12::pow(proof.Cm_ef,ef,r_f);
    pairing(er1,h1,X);
    Fp12::pow(proof.Cm_er1,er1,r_r1);
    pairing(ed1,h1,h_G2);
    Fp12::pow(proof.Cm_ed1,ed1,r_d1);
    pairing(et,h0,h_G2);
    Fp12::pow(proof.Cm_et,et,r_t);
    pairing(eu,h1,h_G2);
    Fp12::pow(proof.Cm_eu,eu,r_u);
    pairing(ev,h2,h_G2);
    Fp12::pow(proof.Cm_ev,ev,r_v);



    //gen Z
    proof.Z_r1=r_r1+c*r1;
    proof.Z_r2=r_r2+c*r2;
    proof.Z_d1=r_d1+c*delta1;
    proof.Z_d2=r_d2+c*delta2;
    proof.Z_t=r_t+c*user_out.t;
    proof.Z_u=r_u+c*user_out.u;
    proof.Z_f=r_f+c*user_out.f;
    proof.Z_v=r_v+c*user_out.v;

    return proof;
}

bool SPK_verify(const SPK_proof &proof, const mcl::bn::G1& h, const mcl::bn::G1 &h0, const mcl::bn::G1& h1,const mcl::bn::G1& h2,const mcl::bn::G2 &h_G2,const mcl::bn::Fp12 &G, const mcl::bn::Fp12 &H,const mcl::bn::Fp12 &H1) {

    G1 lB1,rB1,B1c,h0Z_r1,h1Z_r2;
    G1::mul(B1c,proof.B1,proof.c);
    G1::add(lB1,proof.Cm_B1,B1c);
    G1::mul(h0Z_r1,h0,proof.Z_r1);
    G1::mul(h1Z_r2,h1,proof.Z_r2);
    G1::add(rB1,h0Z_r1,h1Z_r2);
    //std::cout << "[B1 검증] " << (lB1 == rB1 ? "✅ 성공" : "❌ 실패") << std::endl;

    G1 lB1f,rB1f,B1fc,h0Z_d1,h1Z_d2;
    G1::mul(B1fc,proof.B1f,proof.c);
    G1::add(lB1f,proof.Cm_B1f,B1fc);
    G1::mul(rB1f,proof.B1,proof.Z_f);

    //std::cout << "[B1e 검증] " << (lB1f == rB1f ? "✅ 성공" : "❌ 실패") << std::endl;

    Fp12 e_B2X, e_hh_G2, result;

    pairing(e_B2X, proof.B2, X);       // e(A2, w)
    pairing(e_hh_G2, h, h_G2);     // e(g0, h0)

    Fp12 inv_e_hh_G2;
    Fp12::inv(inv_e_hh_G2, e_hh_G2);   // 역원 계산
    Fp12::mul(result, e_B2X, inv_e_hh_G2); // 최종 결과 계산
    
    Fp12 templ1,templ2,templ3,templ4,templ5,el,resultc,inv_Cm_ee;

    Fp12::pow(resultc,result,proof.c);

    Fp12::inv(inv_Cm_ee,proof.Cm_ef);

    Fp12::mul(templ1,inv_Cm_ee,proof.Cm_er1);

    Fp12::mul(templ2,proof.Cm_ed1,templ1);
    Fp12::mul(templ3,proof.Cm_et,templ2);
    Fp12::mul(templ4,proof.Cm_eu,templ3);
    Fp12::mul(templ5,proof.Cm_ev,templ4);
    Fp12::mul(el,templ5,resultc);

    Fp12 tempr1,tempr2,tempr3,tempr4,tempr5,tempr6,inv_tempr1,er;

    Fp12 ef,er1,ed1,et,eu,ev;

    pairing(ef,proof.B2,h_G2);
    pairing(er1,h1,X);
    pairing(ed1,h1,h_G2);
    pairing(et,h0,h_G2);
    pairing(eu,h1,h_G2);
    pairing(ev,h2,h_G2);
    

    Fp12::pow(tempr1,ef,proof.Z_f);
    Fp12::inv(inv_tempr1, tempr1);   // 역원 계산
    
    Fp12::pow(tempr2,er1,proof.Z_r1);
    Fp12::pow(tempr3,ed1,proof.Z_d1);
    Fp12::pow(tempr4,et,proof.Z_t);
    Fp12::pow(tempr5,eu,proof.Z_u);
    Fp12::pow(tempr6,ev,proof.Z_v);

    Fp12 t1,t2,t3,t4;
    Fp12::mul(t1,inv_tempr1,tempr2);
    Fp12::mul(t2,t1,tempr3);
    Fp12::mul(t3,t2,tempr4);
    Fp12::mul(t4,t3,tempr5);
    Fp12::mul(er,t4,tempr6);

    //std::cout << "[페어링 검증] " << (er == el ? "✅ 성공" : "❌ 실패") << std::endl;

    //S 검증

    Fp12 S_c,lS,rS;
    Fp12::pow(S_c,proof.S,proof.c);
    Fp12::mul(lS, S_c, proof.Cm_S);

    Fp12::pow(rS,H ,proof.Z_v);

    //std::cout << "[S 검증] " << (lS==rS ? "✅ 성공" : "❌ 실패") << std::endl;

    //D 검증

    Fp12 D_c,lD,rD,temp_rD1,temp_rD2;
    Fp12::pow(D_c,proof.D,proof.c);
    Fp12::mul(lD, D_c, proof.Cm_D);

    Fp12::pow(temp_rD1,G ,proof.Z_u);
    Fp12::pow(temp_rD2,H1 ,proof.R);
    Fp12::pow(temp_rD2,temp_rD2 ,proof.Z_v);
    Fp12::mul(rD,temp_rD1,temp_rD2);

    //std::cout << "[D 검증] " << (lD==rD ? "✅ 성공" : "❌ 실패") << std::endl;



    return (lB1 == rB1)&&(lB1f == rB1f)&&(er == el)&&(lS==rS)&&(lD==rD); 
    }
bool Payment(User& payer, User& payee,
    const mcl::bn::G1& g, const mcl::bn::G1& g0, const mcl::bn::G1& g1,
    const mcl::bn::G2& g_G2, const mcl::bn::Fp12& G,
    const mcl::bn::G1& h, const mcl::bn::G1& h0, const mcl::bn::G1& h1,
    const mcl::bn::G1& h2, const mcl::bn::G2& h_G2,
    const mcl::bn::Fp12& H, const mcl::bn::Fp12& H1,const std::string INFO){
    //payee side 1 start
    // Step 1: pk_3 proof 생성
    pk_3_proof pk3 = pk_3_prove(payee, g, g0, g1, g_G2, G, INFO);

    // Step 2: pk_3 검증
    if (!pk_3_verify(pk3, g, g0, g1, g_G2, G)) {
        std::cout << "❌ Verification failed: pk_3 is not valid" << std::endl;
        return false;
    }
    //payee side 1 end

    //payer side 1 start
    // Step 3: SPK proof 생성
    SPK_proof spk = SPK_prove(payer, h, h0, h1, h2, h_G2, G, H, H1, INFO, pk3);
    //payer side 1 end

    // Step 4: SPK 검증
    //payee side 2 start
    if (!SPK_verify(spk, h, h0, h1, h2, h_G2, G, H, H1)) {
        std::cout << "❌ Verification failed: SPK is not valid" << std::endl;
        return false;
    }

    // Step 5: 증명 통과 → 저장 및 성공 메시지
    payee.spk = spk;
    //std::cout << "✅ Payment completed successfully." << std::endl;
    return true;
    //payee side 2 end
    }


// Prover: 증명 생성
pk_4_proof pk_4_prove(User &user_out,const mcl::bn::Fr &tp,const mcl::bn::Fr &vp,
    const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2, const mcl::bn::G1 &Cmp,const mcl::bn::Fp12 &G) {

    pk_4_proof proof;

    proof.N=user_out.N;
    proof.M=user_out.M;

    Fr r_tp, r_u,r_vp;
    r_tp.setByCSPRNG();
    r_u.setByCSPRNG();
    r_vp.setByCSPRNG();

    std::string input = proof.INFO + proof.N.getStr(10);
    Fp12 G_t_info;
    hashToGT(G_t_info,input);

    G1::mul(proof.Cm_tp, h0, r_tp);
    G1::mul(proof.Cm_u, h1, r_u);
    G1::mul(proof.Cm_vp,h2,r_vp);

    Fp12::pow(proof.M,G_t_info,user_out.u);
    
    Fp12::pow(proof.Cm_M,G_t_info,r_u);

    proof.c.setByCSPRNG();
    proof.Z_tp = r_tp + proof.c * tp;
    proof.Z_u  = r_u  + proof.c * user_out.u;
    proof.Z_vp = r_vp + proof.c * vp;
    

    return proof;
}

// Verifier: 증명 검증
bool pk_4_verify(const pk_4_proof &proof, const mcl::bn::G1 &Cmp ,const mcl::bn::G1 &h0, const mcl::bn::G1 &h1 , const mcl::bn::G1 &h2 ,const mcl::bn::Fp12 &G) {

    G1 lhs_h0, lhs_h1, lhs_h2,temp, lhs;
    G1::mul(lhs_h0, h0, proof.Z_tp);
    G1::mul(lhs_h1, h1, proof.Z_u);
    G1::mul(lhs_h2, h2, proof.Z_vp);
    G1::add(temp,lhs_h0,lhs_h1);
    G1::add(lhs,temp,lhs_h2);

    G1 rhs_Cmp, rhs_temp1, rhs_temp2,rhs;
    G1::mul(rhs_Cmp, Cmp, proof.c);
    G1::add(rhs_temp1, proof.Cm_tp, proof.Cm_u);
    G1::add(rhs_temp2, rhs_temp1, rhs_Cmp);
    G1::add(rhs,rhs_temp2,proof.Cm_vp);
    //std::cout << "[Cm' 검증] " << (lhs==rhs ? "✅ 성공" : "❌ 실패") << std::endl;

    std::string input = proof.INFO + proof.N.getStr(10);
    Fp12 G_t_info;
    hashToGT(G_t_info,input);

    Fp12 M_c,lM;
    Fp12::pow(M_c,proof.M,proof.c);
    Fp12::mul(lM, M_c, proof.Cm_M);

    Fp12 rM;
    Fp12::pow(rM,G_t_info ,proof.Z_u);

    //std::cout << "[M 검증] " << (lM==rM ? "✅ 성공" : "❌ 실패") << std::endl;
    
    return (lhs == rhs)&&(lM==rM ) ;
}
bool randomise(User &user,
    const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1, const mcl::bn::G2 &g_G2, const mcl::bn::Fp12 &G,
    const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2, const mcl::bn::G2 &h_G2, const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1)
{
    //payee side 1 start
    // Step 0: 검증 전 증명 생성
    pk_3_proof pk3 = pk_3_prove(user, g, g0, g1, g_G2, G,user.INFO);
    if (!pk_3_verify(pk3, g, g0, g1, g_G2, G)) {
        std::cout << "❌ Verification failed: Pk_4 is not valid" << std::endl;
        return false;
    }
    //payee side 1 end

    //bank side 1 start

    if (!SPK_verify(user.spk, h, h0, h1, h2, h_G2, G, H, H1)) {
        std::cout << "❌ Verification failed: SPK is not valid" << std::endl;
        return false;
    }

    // Step 0.5: Double Spending Check
    //find cheater
    /*if (checkDoubleSpending(user.spk)) {//DB에 가서 spk의 s가 존재하는지 안 존재하는지 판단단
        std::cout << "❌ Coin rejected: Double spending detected." << std::endl;
        return false;
    }*/
    

    // Step 0.6: Double Spending 아니라면 DB에 등록
    insertSPKProof(db, user.spk);
    //bank side 1 end

    //payee side 2 start
    // Step 1: User generates commitment Cm
    Fr tp, vp;
    tp.setByCSPRNG();
    vp.setByCSPRNG();

    G1 h0_tp, h1_u, h2_vp, temp, Cmp;
    G1::mul(h0_tp, h0, tp);
    G1::mul(h1_u, h1, user.u);
    G1::mul(h2_vp, h2, vp);
    G1::add(temp, h0_tp, h1_u);
    G1::add(Cmp, h2_vp, temp);

    // Step 2: Prove knowledge of Cm'
    pk_4_proof pk4 = pk_4_prove(user,tp,  vp, h0, h1, h2, Cmp, G);
    //payee side 2 end

    //bank side 2 start
    if (!pk_4_verify(pk4, Cmp, h0, h1, h2, G)) {
        std::cout << "❌ Verification failed: Pk_2 is not valid" << std::endl;
        return false;
    }

    // Step 3: Bank signs the new coin
    Fr tpp, vpp, f;
    tpp.setByCSPRNG();
    vpp.setByCSPRNG();
    f.setByCSPRNG();

    G1 h0_tpp, h2_vpp;
    G1::mul(h0_tpp, h0, tpp);
    G1::mul(h2_vpp, h2, vpp);
    G1::add(temp, h, h0_tpp);
    G1::add(temp, h2_vpp, temp);
    G1::add(temp, Cmp, temp);

    Fr denom = f + beta;
    Fr inv_denom;
    Fr::inv(inv_denom, denom);
    G1::mul(user.B, temp, inv_denom);
    //bank side 2 end

    //payee side 3 start

    // Step 4: Update user secret
    user.f = f;
    user.t = tp + tpp;
    user.v = vp + vpp;

    // Step 5: Verify randomized coin
    G1 h0t, h2v, h0th1u, h0th1uh2v, hh0th1uh2v;
    G1::mul(h0t, h0, user.t);
    G1::mul(h2v, h2, user.v);
    G1::add(h0th1u, h0t, h1_u);
    G1::add(h0th1uh2v, h0th1u, h2v);
    G1::add(hh0th1uh2v, h0th1uh2v, h);

    G2 X, hf, Xhf;
    G2::mul(X, h_G2, beta);
    G2::mul(hf, h_G2, f);
    G2::add(Xhf, hf, X);

    Fp12 e1, e2;
    pairing(e1, user.B, Xhf);
    pairing(e2, hh0th1uh2v, h_G2);

    if (e1 == e2) {
        std::cout << "✅ Coin randomise successful!" << std::endl;
        return true;
    } else {
        std::cout << "❌ Coin randomise failed!" << std::endl;
        return false;
    }
    //payee side 3 end
}
pk_5_proof pk_5_prove(User &user_out,const mcl::bn::Fp12 &G) {

    pk_5_proof proof;

    proof.N=user_out.N;
    proof.M=user_out.M;
    proof.U=user_out.U;

    Fr  r_u;
   
    r_u.setByCSPRNG();
   
    Fp12::pow(proof.Cm_U,G,r_u);

    std::string input = proof.INFO + proof.N.getStr(10);
    Fp12 G_t_info;
    hashToGT(G_t_info,input);

    Fp12::pow(proof.M,G_t_info,user_out.u);
    
    Fp12::pow(proof.Cm_M,G_t_info,r_u);

    proof.c.setByCSPRNG();
  
    proof.Z_u  = r_u  + proof.c * user_out.u;
    

    return proof;
}

// Verifier: 증명 검증
bool pk_5_verify(const pk_5_proof &proof,const mcl::bn::Fp12 &G) {
   
    Fp12 U_c,lU;
    Fp12::pow(U_c,proof.U,proof.c);
    Fp12::mul(lU, U_c, proof.Cm_U);
    Fp12 rU;
    Fp12::pow(rU,G,proof.Z_u);

    //std::cout << "[U 검증] " << (lU==rU ? "✅ 성공" : "❌ 실패") << std::endl;

    std::string input = proof.INFO + proof.N.getStr(10);
    Fp12 G_t_info;
    hashToGT(G_t_info,input);

    Fp12 M_c,lM;
    Fp12::pow(M_c,proof.M,proof.c);
    Fp12::mul(lM, M_c, proof.Cm_M);

    Fp12 rM;
    Fp12::pow(rM,G_t_info ,proof.Z_u);

    //std::cout << "[M 검증] " << (lM==rM ? "✅ 성공" : "❌ 실패") << std::endl;
    
    return (lU == rU)&&(lM==rM ) ;
}

bool finalise(User &user,
    const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1,const mcl::bn::G2 &g_G2,const mcl::bn::Fp12 &G,
    const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2,const mcl::bn::G2 &h_G2,const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1) {
// Step 0: 검증 전 증명 생성
//payee side 1 start
pk_5_proof pk5 = pk_5_prove(user,G);
//payee side 1 start

//bank side 1 start
if (!pk_5_verify(pk5, G)) {
std::cout << "❌ Verification failed: Pk_5 is not valid" << std::endl;
return false;}


if (!SPK_verify(user.spk, h, h0, h1, h2, h_G2, G, H, H1)) {
std::cout << "❌ Verification failed: SPK is not valid" << std::endl;
return false;}
//bank side 1 end

}
bool initDatabase(sqlite3* &db) {
    int rc = sqlite3_open("spk_proof.db", &db);
    if (rc) {
        std::cerr << "❌ Cannot open DB: " << sqlite3_errmsg(db) << std::endl;
        return false;
    }

    const char* sql_create = R"(
        CREATE TABLE IF NOT EXISTS spk_proof (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            B1 TEXT, B1f TEXT, B2 TEXT,
            S TEXT, D TEXT,
            Cm_B1 TEXT, Cm_B1f TEXT,
            Cm_ef TEXT, Cm_er1 TEXT, Cm_ed1 TEXT,
            Cm_et TEXT, Cm_eu TEXT, Cm_ev TEXT,
            Cm_S TEXT, Cm_D TEXT,
            c TEXT, Z_r1 TEXT, Z_r2 TEXT,
            Z_d1 TEXT, Z_d2 TEXT, Z_t TEXT,
            Z_u TEXT, Z_f TEXT, Z_v TEXT,
            R TEXT, N TEXT, M TEXT
        );
    )";

    char* errMsg = nullptr;
    rc = sqlite3_exec(db, sql_create, nullptr, nullptr, &errMsg);
    if (rc != SQLITE_OK) {
        std::cerr << "❌ Table creation failed: " << errMsg << std::endl;
        sqlite3_free(errMsg);
        return false;
    }

    std::cout << "✅ Database initialized.\n";
    return true;
}
//spk를  DB에 저장하는 코드드
bool insertSPKProof(sqlite3* db, const SPK_proof& proof) {
    std::stringstream ss;
    ss << "INSERT INTO spk_proof (B1, B1f, B2, S, D, Cm_B1, Cm_B1f, Cm_ef, Cm_er1, Cm_ed1, "
       << "Cm_et, Cm_eu, Cm_ev, Cm_S, Cm_D, c, Z_r1, Z_r2, Z_d1, Z_d2, Z_t, Z_u, Z_f, Z_v, R, N, M) VALUES ('"
       << proof.B1.getStr(16) << "', '"
       << proof.B1f.getStr(16) << "', '"
       << proof.B2.getStr(16) << "', '"
       << proof.S.getStr(16) << "', '"
       << proof.D.getStr(16) << "', '"
       << proof.Cm_B1.getStr(16) << "', '"
       << proof.Cm_B1f.getStr(16) << "', '"
       << proof.Cm_ef.getStr(16) << "', '"
       << proof.Cm_er1.getStr(16) << "', '"
       << proof.Cm_ed1.getStr(16) << "', '"
       << proof.Cm_et.getStr(16) << "', '"
       << proof.Cm_eu.getStr(16) << "', '"
       << proof.Cm_ev.getStr(16) << "', '"
       << proof.Cm_S.getStr(16) << "', '"
       << proof.Cm_D.getStr(16) << "', '"
       << proof.c.getStr(16) << "', '"
       << proof.Z_r1.getStr(16) << "', '"
       << proof.Z_r2.getStr(16) << "', '"
       << proof.Z_d1.getStr(16) << "', '"
       << proof.Z_d2.getStr(16) << "', '"
       << proof.Z_t.getStr(16) << "', '"
       << proof.Z_u.getStr(16) << "', '"
       << proof.Z_f.getStr(16) << "', '"
       << proof.Z_v.getStr(16) << "', '"
       << proof.R.getStr(16) << "', '"
       << proof.N.getStr(16) << "', '"
       << proof.M.getStr(16) << "');";

    std::string sql = ss.str();
    char* errMsg = nullptr;
    int rc = sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &errMsg);
    if (rc != SQLITE_OK) {
        std::cerr << "❌ Insert failed: " << errMsg << std::endl;
        sqlite3_free(errMsg);
        return false;
    }

    //std::cout << "✅ SPK_proof inserted into DB." << std::endl;
    return true;
}
//DB가서 더블스펜딩이 일어났나 안일어났나 확인하는 코드
bool checkDoubleSpending(const SPK_proof& proof) {
    std::string sStr = proof.S.getStr(16);

    std::string query = "SELECT COUNT(*) FROM spk_proof WHERE S = '" + sStr + "';";
    sqlite3_stmt* stmt;

    int rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        std::cerr << "❌ Failed to prepare statement: " << sqlite3_errmsg(db) << std::endl;
        return false;
    }

    rc = sqlite3_step(stmt);
    if (rc == SQLITE_ROW) {
        int count = sqlite3_column_int(stmt, 0);
        sqlite3_finalize(stmt);

        if (count > 0) {
            std::cout << "🚨 Double Spending Detected! Matching S exists in DB." << std::endl;
            return true;
        } else {
            std::cout << "✅ No double spending detected." << std::endl;
            return false;
        }
    } else {
        std::cerr << "❌ Failed to execute statement." << std::endl;
        sqlite3_finalize(stmt);
        return false;
    }
}

bool account_est_time(User &user_out, mcl::bn::G1& g, mcl::bn::G1& g0, mcl::bn::G1& g1, mcl::bn::G2& g_G2,const mcl::bn::Fp12& G,
    long long& user1_time, long long& bank_time, long long& user2_time) {

auto t1 = high_resolution_clock::now();

// --- User side 1 ---
Fr sp;
sp.setByCSPRNG();
G1 g0sp, g1u, Cm;
G1::mul(g0sp, g0, sp);
G1::mul(g1u, g1, user_out.u);
G1::add(Cm, g0sp, g1u);

pk_1_proof proof = pk_1_prove(sp, user_out.u, g0, g1, Cm, G);

auto t2 = high_resolution_clock::now();
// --- Bank side ---
if (!pk_1_verify(proof, Cm, g0, g1, G)) {
std::cout << "❌ Verification failed: Pk_1 is not valid" << std::endl;
return false;
}




Fr spp, e;
spp.setByCSPRNG();
e.setByCSPRNG();

G1 g0spp, gg0spp, Cm_gg0spp;
G1::mul(g0spp, g0, spp);
G1::add(gg0spp, g, g0spp);
G1::add(Cm_gg0spp, Cm, gg0spp);

Fr denom = e + alpha;
Fr inv_denom;
Fr::inv(inv_denom, denom);

G1 A;
G1::mul(A, Cm_gg0spp, inv_denom);

auto t3 = high_resolution_clock::now();

// --- User side 2 ---
user_out.A = A;
user_out.e = e;
user_out.s = sp + spp;

G2 W, he, Whe;
G2::mul(W, g_G2, alpha);
G2::mul(he, g_G2, e);
G2::add(Whe, he, W);

G1 g0s, gg0s, gg0sg1u;
G1::mul(g0s, g0, user_out.s);
G1::add(gg0s, g, g0s);
G1::add(gg0sg1u, gg0s, g1u);

Fp12 e1, e2;
pairing(e1, user_out.A, Whe);
pairing(e2, gg0sg1u, g_G2);

auto t4 = high_resolution_clock::now();

// Save timings
user1_time = duration_cast<nanoseconds>(t2 - t1).count();
bank_time  = duration_cast<nanoseconds>(t3 - t2).count();
user2_time = duration_cast<nanoseconds>(t4 - t3).count();

return e1 == e2;
}
bool withdraw_coin_time(User &user_out, mcl::bn::G1& h, mcl::bn::G1& h0, mcl::bn::G1& h1, mcl::bn::G1& h2,
    mcl::bn::G2& h_G2, const mcl::bn::Fp12& G,
    long long& user1_time, long long& bank_time, long long& user2_time) {

auto t1 = high_resolution_clock::now();

// --- User side 1 ---
Fr tp, vp;
tp.setByCSPRNG();
vp.setByCSPRNG();

G1 h0_tp, h1_u, h2_vp, temp, Cmp;
G1::mul(h0_tp, h0, tp);
G1::mul(h1_u, h1, user_out.u);
G1::mul(h2_vp, h2, vp);
G1::add(temp, h0_tp, h1_u);
G1::add(Cmp, h2_vp, temp);

pk_2_proof proof = pk_2_prove(tp, user_out.u, vp, h0, h1, h2, Cmp, G);

auto t2 = high_resolution_clock::now();
// --- Bank side ---
if (!pk_2_verify(proof, Cmp, h0, h1, h2, G)) {
std::cout << "❌ Verification failed: Pk_2 is not valid" << std::endl;
return false;
}




G1 h0_tpp, h2_vpp;
Fr tpp, vpp, f;
tpp.setByCSPRNG();
vpp.setByCSPRNG();
f.setByCSPRNG();

G1::mul(h0_tpp, h0, tpp);
G1::mul(h2_vpp, h2, vpp);
G1::add(temp, h, h0_tpp);
G1::add(temp, h2_vpp, temp);
G1::add(temp, Cmp, temp);

Fr denom = f + beta;
Fr inv_denom;
Fr::inv(inv_denom, denom);
G1::mul(user_out.B, temp, inv_denom);

auto t3 = high_resolution_clock::now();

// --- User side 2 ---
user_out.f = f;
user_out.t = tp + tpp;
user_out.v = vp + vpp;

G2 X, hf, Xhf;
G2::mul(X, h_G2, beta);
G2::mul(hf, h_G2, f);
G2::add(Xhf, hf, X);

G1 h0t, h2v, h0th1u, h0th1uh2v, hh0th1uh2v;
G1::mul(h0t, h0, user_out.t);
G1::mul(h2v, h2, user_out.v);
G1::add(h0th1u, h0t, h1_u);
G1::add(h0th1uh2v, h0th1u, h2v);
G1::add(hh0th1uh2v, h0th1uh2v, h);

Fp12 e1, e2;
pairing(e1, user_out.B, Xhf);
pairing(e2, hh0th1uh2v, h_G2);

auto t4 = high_resolution_clock::now();

// 시간 기록
user1_time = duration_cast<nanoseconds>(t2 - t1).count();
bank_time  = duration_cast<nanoseconds>(t3 - t2).count();
user2_time = duration_cast<nanoseconds>(t4 - t3).count();

return (e1 == e2);
}

bool Payment_time(User& payer, User& payee,
    const mcl::bn::G1& g, const mcl::bn::G1& g0, const mcl::bn::G1& g1,
    const mcl::bn::G2& g_G2, const mcl::bn::Fp12& G,
    const mcl::bn::G1& h, const mcl::bn::G1& h0, const mcl::bn::G1& h1,
    const mcl::bn::G1& h2, const mcl::bn::G2& h_G2,
    const mcl::bn::Fp12& H, const mcl::bn::Fp12& H1,
    const std::string INFO,
    long long& payee1_time, long long& payer_time, long long& payee2_time) 
{
    auto t1 = high_resolution_clock::now();

    // --- Payee side 1 ---
    pk_3_proof pk3 = pk_3_prove(payee, g, g0, g1, g_G2, G, INFO);

    auto t2 = high_resolution_clock::now();
    // --- Payer side 1 ---
    if (!pk_3_verify(pk3, g, g0, g1, g_G2, G)) {
        std::cout << "❌ Verification failed: pk_3 is not valid" << std::endl;
        return false;
    }
    
    SPK_proof spk = SPK_prove(payer, h, h0, h1, h2, h_G2, G, H, H1, INFO, pk3);

    auto t3 = high_resolution_clock::now();

    // --- Payee side 2 ---
    if (!SPK_verify(spk, h, h0, h1, h2, h_G2, G, H, H1)) {
        std::cout << "❌ Verification failed: SPK is not valid" << std::endl;
        return false;
    }

    payee.spk = spk;

    auto t4 = high_resolution_clock::now();

    payee1_time = duration_cast<nanoseconds>(t2 - t1).count();
    payer_time  = duration_cast<nanoseconds>(t3 - t2).count();
    payee2_time = duration_cast<nanoseconds>(t4 - t3).count();

    
}
bool randomise_time(User &user,
    const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1, const mcl::bn::G2 &g_G2, const mcl::bn::Fp12 &G,
    const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2, const mcl::bn::G2 &h_G2,
    const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1,
    long long& payee1_time, long long& bank1_time, long long& payee2_time, long long& bank2_time, long long& payee3_time)
{
    auto t1 = high_resolution_clock::now();

    // --- Payee Side 1 ---
    pk_3_proof pk3 = pk_3_prove(user, g, g0, g1, g_G2, G, user.INFO);

    auto t2 = high_resolution_clock::now();
    if (!pk_3_verify(pk3, g, g0, g1, g_G2, G)) {
        std::cout << "❌ Verification failed: Pk_4 is not valid" << std::endl;
        return false;
    }

    

    // --- Bank Side 1 ---
    if (!SPK_verify(user.spk, h, h0, h1, h2, h_G2, G, H, H1)) {
        std::cout << "❌ Verification failed: SPK is not valid" << std::endl;
        return false;
    }

    insertSPKProof(db, user.spk);  // Double spending check 및 등록

    auto t3 = high_resolution_clock::now();

    // --- Payee Side 2 ---
    Fr tp, vp;
    tp.setByCSPRNG();
    vp.setByCSPRNG();

    G1 h0_tp, h1_u, h2_vp, temp, Cmp;
    G1::mul(h0_tp, h0, tp);
    G1::mul(h1_u, h1, user.u);
    G1::mul(h2_vp, h2, vp);
    G1::add(temp, h0_tp, h1_u);
    G1::add(Cmp, h2_vp, temp);

    pk_4_proof pk4 = pk_4_prove(user, tp, vp, h0, h1, h2, Cmp, G);

    auto t4 = high_resolution_clock::now();

    // --- Bank Side 2 ---
    if (!pk_4_verify(pk4, Cmp, h0, h1, h2, G)) {
        std::cout << "❌ Verification failed: Pk_2 is not valid" << std::endl;
        return false;
    }

    Fr tpp, vpp, f;
    tpp.setByCSPRNG();
    vpp.setByCSPRNG();
    f.setByCSPRNG();

    G1 h0_tpp, h2_vpp;
    G1::mul(h0_tpp, h0, tpp);
    G1::mul(h2_vpp, h2, vpp);
    G1::add(temp, h, h0_tpp);
    G1::add(temp, h2_vpp, temp);
    G1::add(temp, Cmp, temp);

    Fr denom = f + beta;
    Fr inv_denom;
    Fr::inv(inv_denom, denom);
    G1::mul(user.B, temp, inv_denom);

    auto t5 = high_resolution_clock::now();

    // --- Payee Side 3 ---
    user.f = f;
    user.t = tp + tpp;
    user.v = vp + vpp;

    G1 h0t, h2v, h0th1u, h0th1uh2v, hh0th1uh2v;
    G1::mul(h0t, h0, user.t);
    G1::mul(h2v, h2, user.v);
    G1::add(h0th1u, h0t, h1_u);
    G1::add(h0th1uh2v, h0th1u, h2v);
    G1::add(hh0th1uh2v, h0th1uh2v, h);

    G2 X, hf, Xhf;
    G2::mul(X, h_G2, beta);
    G2::mul(hf, h_G2, f);
    G2::add(Xhf, hf, X);

    Fp12 e1, e2;
    pairing(e1, user.B, Xhf);
    pairing(e2, hh0th1uh2v, h_G2);

    if (e1 == e2) {
        //std::cout << "✅ Coin randomise successful!" << std::endl;
        //return true;
    } else {
        //std::cout << "❌ Coin randomise failed!" << std::endl;
        return false;
    }

    auto t6 = high_resolution_clock::now();

    payee1_time = duration_cast<nanoseconds>(t2 - t1).count();
    bank1_time  = duration_cast<nanoseconds>(t3 - t2).count();
    payee2_time = duration_cast<nanoseconds>(t4 - t3).count();
    bank2_time  = duration_cast<nanoseconds>(t5 - t4).count();
    payee3_time = duration_cast<nanoseconds>(t6 - t5).count();
    return true;
    }
bool finalise_time(User &user,
    const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1, const mcl::bn::G2 &g_G2, const mcl::bn::Fp12 &G,
    const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2, const mcl::bn::G2 &h_G2,
    const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1,
    long long &payee1_time, long long &bank1_time)
{
    auto t1 = high_resolution_clock::now();

    // --- Payee side 1 ---
    pk_5_proof pk5 = pk_5_prove(user, G);

    auto t2 = high_resolution_clock::now();

    // --- Bank side 1 ---
    if (!pk_5_verify(pk5, G)) {
        std::cout << "❌ Verification failed: Pk_5 is not valid" << std::endl;
        return false;
    }

    if (!SPK_verify(user.spk, h, h0, h1, h2, h_G2, G, H, H1)) {
        std::cout << "❌ Verification failed: SPK is not valid" << std::endl;
        return false;
    }

    auto t3 = high_resolution_clock::now();

    payee1_time = duration_cast<nanoseconds>(t2 - t1).count();
    bank1_time = duration_cast<nanoseconds>(t3 - t2).count();

}