#pragma once
#include <mcl/bn.hpp>  // 반드시 이 라인이 필요합니다
using namespace mcl::bn;
#include <sqlite3.h>
#include <chrono>  // 이 헤더 추가

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdint.h>
#include <cstring>

extern mcl::bn::Fr alpha;
extern mcl::bn::Fr beta;
extern mcl::bn::Fr x_GE;
extern mcl::bn::Fr x_h;
extern mcl::bn::G2 W;
extern mcl::bn::G2 X;
extern mcl::bn::Fp12 P;
extern mcl::bn::G1 Q;
extern sqlite3* db;  // 다른 파일에서 선언된 전역 DB 포인터 사용
// 사용자 정보를 저장하는 구조체



// 시그마 프로토콜 증명 구조체
struct pk_1_proof {

    mcl::bn::Fp12 U;
    mcl::bn::G1 Cm_sp;
    mcl::bn::G1 Cm_u;
    mcl::bn::Fp12 Cm_U;
    mcl::bn::Fr c;
    mcl::bn::Fr Z_sp;
    mcl::bn::Fr Z_u;
   
};
struct pk_2_proof {
    mcl::bn::Fp12 U;
    mcl::bn::G1 Cm_tp;
    mcl::bn::G1 Cm_u;
    mcl::bn::G1 Cm_vp;
    mcl::bn::Fp12 Cm_U;
    mcl::bn::Fr c;//
    mcl::bn::Fr Z_tp;
    mcl::bn::Fr Z_u;
    mcl::bn::Fr Z_vp;
    mcl::bn::Fr Z_U;
};
struct pk_3_proof {
    mcl::bn::Fr N;
    mcl::bn::Fp12 M;
    std::string INFO;
    
    mcl::bn::G1 A1;
    mcl::bn::G1 A1e;
    mcl::bn::G1 A2;
    
    mcl::bn::G1 Cm_A1;
    mcl::bn::G1 Cm_A1e;
    mcl::bn::Fp12 Cm_ee;//e(A2,g)^r_e
    mcl::bn::Fp12 Cm_er1;//e(g1,g)^r_r1
    mcl::bn::Fp12 Cm_ed1;//e(g1,g)^r_d1
    mcl::bn::Fp12 Cm_es;//e(g0,g)^r_s
    mcl::bn::Fp12 Cm_eu;//e(g1,g)^r_u
    mcl::bn::Fp12 Cm_M;
    mcl::bn::Fr c;//
    mcl::bn::Fr Z_r1;
    mcl::bn::Fr Z_r2;
    mcl::bn::Fr Z_d1;
    mcl::bn::Fr Z_d2;
    mcl::bn::Fr Z_s;
    mcl::bn::Fr Z_u;
    mcl::bn::Fr Z_e;
};

struct pk_4_proof {
    
    mcl::bn::Fr N;
    mcl::bn::Fp12 M;
    
    mcl::bn::Fp12 U;
    std::string INFO;

    mcl::bn::G1 Cm_tp;
    mcl::bn::G1 Cm_u;
    mcl::bn::G1 Cm_vp;
    mcl::bn::Fp12 Cm_M;
    mcl::bn::Fr c;//
    mcl::bn::Fr Z_tp;
    mcl::bn::Fr Z_u;
    mcl::bn::Fr Z_vp;
    mcl::bn::Fr Z_M;
};
struct SPK1_proof {
    mcl::bn::G1 B1;
    mcl::bn::G1 B1f;
    mcl::bn::G1 B2;

    mcl::bn::G1 Cm_B1;
    mcl::bn::G1 Cm_B1f;
    mcl::bn::Fp12 Cm_Bef;//e(B2,h)^r_e
    mcl::bn::Fp12 Cm_Ber1;//e(h1,h)^r_r1
    mcl::bn::Fp12 Cm_Beg1;//e(h1,h)^r_g1
    mcl::bn::Fp12 Cm_Bet;//e(h0,h)^r_s
    mcl::bn::Fp12 Cm_Beu;//e(h1,h)^r_u
    mcl::bn::Fp12 Cm_Bev;//e(h0,h)^r_s

    //Tracing tag
    mcl::bn::G1 T;
    mcl::bn::G1 Cm_T;
    //encrypted backward tag
    mcl::bn::G1 backward_C1;
    mcl::bn::G1 backward_C2;
    mcl::bn::G1 Cm_backward_r1;
    mcl::bn::G1 Cm_backward_C1;
    mcl::bn::G1 Cm_backward_C2;

    mcl::bn::Fp12 S;
    mcl::bn::Fp12 D;

    mcl::bn::Fp12 Cm_S;
    mcl::bn::Fp12 Cm_D;
    //for user identifier
    mcl::bn::Fp12 userID_C1;
    mcl::bn::Fp12 userID_C2;

    mcl::bn::Fp12 Cm_GE;//for userID_C1
    mcl::bn::Fp12 Cm_G;//for userID_C2
    mcl::bn::Fp12 Cm_P;//for userID_C2

    mcl::bn::Fr c;//
    mcl::bn::Fr Z_r1;
    mcl::bn::Fr Z_r2;
    mcl::bn::Fr Z_g1;
    mcl::bn::Fr Z_g2;
    mcl::bn::Fr Z_t;
    mcl::bn::Fr Z_f;
    mcl::bn::Fr Z_v;

    mcl::bn::Fr Z_d1;
    mcl::bn::Fr Z_d2;
    mcl::bn::Fr Z_s;
    mcl::bn::Fr Z_u;
    mcl::bn::Fr Z_e;
    //for gen R
    mcl::bn::Fr R;
    mcl::bn::Fr N;
    mcl::bn::Fp12 M;
};
struct SPK2_proof {
    mcl::bn::Fr N;
    mcl::bn::Fp12 M;
    std::string INFO;
    
    mcl::bn::G1 A1;
    mcl::bn::G1 A1e;
    mcl::bn::G1 A2;
    
    mcl::bn::G1 Cm_A1;
    mcl::bn::G1 Cm_A1e;
    mcl::bn::Fp12 Cm_ee;//e(A2,g)^r_e
    mcl::bn::Fp12 Cm_er1;//e(g1,g)^r_r1
    mcl::bn::Fp12 Cm_ed1;//e(g1,g)^r_d1
    mcl::bn::Fp12 Cm_es;//e(g0,g)^r_s
    mcl::bn::Fp12 Cm_eu;//e(g1,g)^r_u
    mcl::bn::Fp12 Cm_M;

    //for user identifier
    mcl::bn::Fp12 userID_C1;
    mcl::bn::Fp12 userID_C2;

    mcl::bn::Fp12 Cm_GE;//for userID_C1
    mcl::bn::Fp12 Cm_G;//for userID_C2
    mcl::bn::Fp12 Cm_P;//for userID_C2

    mcl::bn::Fr c;//
    mcl::bn::Fr Z_r1;
    mcl::bn::Fr Z_r2;
    mcl::bn::Fr Z_d1;
    mcl::bn::Fr Z_d2;
    mcl::bn::Fr Z_s;
    mcl::bn::Fr Z_u;
    mcl::bn::Fr Z_e;
};
struct SPK3_proof {
    

    mcl::bn::Fr c;//
    mcl::bn::Fr Z_r1;
    mcl::bn::Fr Z_t;
    mcl::bn::Fr Z_u;
    mcl::bn::Fr Z_v;

    mcl::bn::G1 Cm_h0;
    mcl::bn::G1 Cm_h1;
    mcl::bn::G1 Cm_h2;

    //for gen R
    mcl::bn::G1 forward_C1;
    mcl::bn::G1 forward_C2;
    mcl::bn::G1 Cm_forward_r1;
    mcl::bn::G1 Cm_forward_C1;
    mcl::bn::G1 Cm_forward_C2;
};
struct User {
    mcl::bn::G1 A;     // 서명
    mcl::bn::Fr e;     // 서명의 일부
    mcl::bn::Fr s;     // sp + spp
    mcl::bn::Fr u;     // 랜덤값 (계정 생성용)
    mcl::bn::G1 B;     // 서명
    mcl::bn::Fr f;     // 서명의 일부
    mcl::bn::Fr t;     // sp + spp
    mcl::bn::Fr v;   
    mcl::bn::Fp12 U;   
    mcl::bn::Fr N;
    mcl::bn::Fp12 M;
    std::string INFO;
    SPK1_proof spk1;     // 사용자 관련 SPK 증명 정보
    SPK2_proof spk2;
    SPK3_proof spk3;
};
// 함수 선언
// 함수 선언
pk_1_proof pk_1_prove(const Fr &sp, const Fr &u, const G1 &g0, const G1 &g1, const G1 &Cm,const Fp12 &G) ;
pk_2_proof pk_2_prove(const mcl::bn::Fr &tp, const mcl::bn::Fr &u,const mcl::bn::Fr &vp,
    const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2, const mcl::bn::G1 &Cmp,const mcl::bn::Fp12 &G);
bool pk_1_verify(const pk_1_proof &proof, const G1 &Cm, const G1 &g0, const G1 &g1,const Fp12 &G) ;
bool pk_2_verify(const pk_2_proof &proof, const mcl::bn::G1 &Cmp, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1 , const mcl::bn::G1 &h2 ,const mcl::bn::Fp12 &G) ;
void setupG_Generators(mcl::bn::G1& g, mcl::bn::G1& g0, mcl::bn::G1& g1, mcl::bn::G1& ge, mcl::bn::G2& gk);
void setupH_Generators(mcl::bn::G1& h, mcl::bn::G1& h0, mcl::bn::G1& h1, mcl::bn::G1& h2,mcl::bn::G1& he, mcl::bn::G1& ht, mcl::bn::G2& g);
bool account_est(User &user_out,mcl::bn::G1& g, mcl::bn::G1& g0, mcl::bn::G1& g1, mcl::bn::G2& g_G2,const mcl::bn::Fp12 &G); 
bool withdraw_coin(User &user_out, mcl::bn::G1& h, mcl::bn::G1& h0, mcl::bn::G1& h1,mcl::bn::G1& h2, mcl::bn::G2& h_G2,const mcl::bn::Fp12& G);
void hashToG1(mcl::bn::G1& g, const std::string& label);
void hashToG2(mcl::bn::G2& h, const std::string& label);
void hashToFr(mcl::bn::Fr& target,const std::string& msg);
void hashToGT(mcl::bn::Fp12& target,const std::string& msg);

pk_3_proof pk_3_prove(User &user_out,const mcl::bn::G1 g, const mcl::bn::G1 &g0, const mcl::bn::G1 g1,
    const mcl::bn::G2 g_G2, const mcl::bn::Fp12 &G,std::string INFO);
bool pk_3_verify(const pk_3_proof &proof, const mcl::bn::G1 g, const mcl::bn::G1 &g0, const mcl::bn::G1 g1,
    const mcl::bn::G2 g_G2, const mcl::bn::Fp12 &G);
SPK1_proof SPK1_prove(User &user_out,const mcl::bn::G1 h, const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2,const mcl::bn::G1 hE,const mcl::bn::G1 ht,const mcl::bn::G2 h_G2,const mcl::bn::Fp12 &H ,const mcl::bn::Fp12 &H1,
    const mcl::bn::Fp12 &G,const mcl::bn::Fp12 &GE,std::string INFO,pk_3_proof pk_3);
bool SPK1_verify(const SPK1_proof &proof, const mcl::bn::G1 h, const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2,const mcl::bn::G1 hE,const mcl::bn::G1 ht,const mcl::bn::G2 h_G2,const mcl::bn::Fp12 &H ,const mcl::bn::Fp12 &H1,
    const mcl::bn::Fp12 &G,const mcl::bn::Fp12 &GE);
bool Payment(User& payer, User& payee,
    const mcl::bn::G1& g, const mcl::bn::G1& g0, const mcl::bn::G1& g1,const mcl::bn::G1& gE,
    const mcl::bn::G2& g_G2, const mcl::bn::Fp12& G,const mcl::bn::Fp12& GE,
    const mcl::bn::G1& h, const mcl::bn::G1& h0, const mcl::bn::G1& h1,
    const mcl::bn::G1& h2,const mcl::bn::G1& hE,const mcl::bn::G1& ht, const mcl::bn::G2& h_G2,
    const mcl::bn::Fp12& H, const mcl::bn::Fp12& H1);
pk_4_proof pk_4_prove(User &user_out,const mcl::bn::Fr &tp,const mcl::bn::Fr &vp,
    const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2, const mcl::bn::G1 &Cmp,const mcl::bn::Fp12 &G);
    bool pk_4_verify(const pk_4_proof &proof, const mcl::bn::G1 &Cmp ,const mcl::bn::G1 &h0, const mcl::bn::G1 &h1 , const mcl::bn::G1 &h2 ,const mcl::bn::Fp12 &G);
SPK2_proof SPK2_prove(User &user_out,const mcl::bn::G1 g, const mcl::bn::G1 &g0, const mcl::bn::G1 g1,const mcl::bn::G2 g_G2, const mcl::bn::Fp12 &G,const mcl::bn::Fp12 &GE) ;
bool SPK2_verify(const SPK2_proof &proof, const mcl::bn::G1 g, const mcl::bn::G1 &g0, const mcl::bn::G1 g1,const mcl::bn::G2 g_G2, const mcl::bn::Fp12 &G,const mcl::bn::Fp12 &GE); 
SPK3_proof SPK3_prove(User &user_out,const mcl::bn::G1 h, const mcl::bn::G1 &h0, const mcl::bn::G1 h1, const mcl::bn::G1 h2,const mcl::bn::G1 hE,const mcl::bn::G1 ht,const mcl::bn::G2 h_G2);
bool SPK3_verify(const SPK3_proof &proof, const mcl::bn::G1& h, const mcl::bn::G1 &h0, const mcl::bn::G1& h1,const mcl::bn::G1& h2,const mcl::bn::G1 hE,const mcl::bn::G1 ht,const mcl::bn::G2 &h_G2,mcl::bn::G1 &bank_B,mcl::bn::Fr &bank_f);
bool randomise(User &user,
    const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1, const mcl::bn::G2 &g_G2, const mcl::bn::Fp12 &G,const mcl::bn::Fp12 &GE,
    const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2,const mcl::bn::G1 &hE,const mcl::bn::G1 &ht, const mcl::bn::G2 &h_G2, const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1);
    bool finalise(User &user,
        const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1, const mcl::bn::G2 &g_G2, const mcl::bn::Fp12 &G,const mcl::bn::Fp12 &GE,
        const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2,const mcl::bn::G1 &hE,const mcl::bn::G1 &ht, const mcl::bn::G2 &h_G2, const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1);
    bool initDatabase(sqlite3*& db) ;
bool insertSPKProof(sqlite3* db,
    const SPK1_proof& spk1,
    const SPK2_proof& spk2,
    const SPK3_proof& spk3,
    const mcl::bn::G1& bank_B) ;
bool insertSPKProofForFinalise(sqlite3* db,
    const SPK1_proof& spk1,
    const SPK2_proof& spk2);
bool checkDoubleSpending(sqlite3* db, const mcl::bn::Fp12& S);
bool insertUserU(sqlite3* db, int user_id, const mcl::bn::Fp12& U) ;
bool account_est_time(User &user_out, mcl::bn::G1& g, mcl::bn::G1& g0, mcl::bn::G1& g1, mcl::bn::G2& g_G2, const mcl::bn::Fp12& G,
    long long& user1_time, long long& bank_time, long long& user2_time);
bool withdraw_coin_time(User &user_out, mcl::bn::G1& h, mcl::bn::G1& h0, mcl::bn::G1& h1, mcl::bn::G1& h2,
    mcl::bn::G2& h_G2, const mcl::bn::Fp12& G,
    long long& user1_time, long long& bank_time, long long& user2_time);
bool Payment_time(User& payer, User& payee,
    const mcl::bn::G1& g, const mcl::bn::G1& g0, const mcl::bn::G1& g1, const mcl::bn::G1& gE,
    const mcl::bn::G2& g_G2, const mcl::bn::Fp12& G, const mcl::bn::Fp12& GE,
    const mcl::bn::G1& h, const mcl::bn::G1& h0, const mcl::bn::G1& h1,
    const mcl::bn::G1& h2, const mcl::bn::G1& hE, const mcl::bn::G1& ht, const mcl::bn::G2& h_G2,
    const mcl::bn::Fp12& H, const mcl::bn::Fp12& H1,
    long long& payee1_time, long long& payer_time, long long& payee2_time);
bool randomise_time(User &user,
    const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1, const mcl::bn::G2 &g_G2, const mcl::bn::Fp12 &G, const mcl::bn::Fp12 &GE,
    const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2, const mcl::bn::G1 &hE, const mcl::bn::G1 &ht, const mcl::bn::G2 &h_G2,
    const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1,
    long long &t1, long long &t2, long long &t3, long long &t4, long long &t5, long long &t6);
bool finalise_time(User &user,
    const mcl::bn::G1 &g, const mcl::bn::G1 &g0, const mcl::bn::G1 &g1, const mcl::bn::G2 &g_G2, const mcl::bn::Fp12 &G, const mcl::bn::Fp12 &GE,
    const mcl::bn::G1 &h, const mcl::bn::G1 &h0, const mcl::bn::G1 &h1, const mcl::bn::G1 &h2, const mcl::bn::G1 &hE, const mcl::bn::G1 &ht,
    const mcl::bn::G2 &h_G2, const mcl::bn::Fp12 &H, const mcl::bn::Fp12 &H1,
    long long &payee_time, long long &bank_time);