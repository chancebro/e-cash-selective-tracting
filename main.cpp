#include <iostream>
#include <sqlite3.h>
#include <mcl/bn.hpp>
#include "tracing.h"
#include <chrono>  // 이 헤더 추가
#include <fstream>

using namespace std::chrono;  // 네임스페이스 생략 위해
using namespace mcl::bn;
void measure_account_est_origin() {
    std::cout << "[account_est] Measurement started...\n";
    if (!initDatabase(db)) return;
    initPairing(); // pairing 초기화

    // Generator 및 pairing 초기화
    G1 g, g0, g1, ge, h, h0, h1, h2, he, ht;
    G2 g_G2, h_G2;
    Fp12 G, H, H1;

    setupH_Generators(h, h0, h1, h2, he, ht, h_G2);
    setupG_Generators(g, g0, g1, ge, g_G2);
    pairing(G, g, g_G2);

    // 사용자 프라이빗 키 초기화용
    User user1, user2;
    user1.u.setByCSPRNG();
    user2.u.setByCSPRNG();
    std::string INFO = "send me a coin";

    Fp12::pow(user1.U, G, user1.u);
    Fp12::pow(user2.U, G, user2.u);

    // 키 생성
    alpha = 1;
    beta  = 2;
    G2::mul(W, g_G2, alpha);
    G2::mul(X, h_G2, beta);

    // 측정 로그 파일 생성
    std::ofstream csv("account_est_timing_origin.csv");
    csv << "Iteration,UserSide1_ns,BankSide_ns,UserSide2_ns\n";

    for (int i = 0; i < 10000; ++i) {
        User user;
        user.u.setByCSPRNG();
        Fp12::pow(user.U, G, user.u);

        long long user1_t = 0, bank_t = 0, user2_t = 0;

        bool success = account_est_time(user, g, g0, g1, g_G2, G,
                                        user1_t, bank_t, user2_t);
        if (!success) {
            std::cerr << "Iteration " << i << ": account_est failed\n";
        }

        csv << i << "," << user1_t << "," << bank_t << "," << user2_t << "\n";
    }

    csv.close();
    sqlite3_close(db);
    std::cout << "[account_est] Measurement completed.\n";
}
void measure_withdraw_coin_origin() {
    std::cout << "[withdraw_coin] Measurement started...\n";
    if (!initDatabase(db)) return;
    initPairing(); // pairing 초기화

    // Generator 및 pairing 초기화
    G1 g, g0, g1, ge, h, h0, h1, h2, he, ht;
    G2 g_G2, h_G2;
    Fp12 G, H, H1;

    setupH_Generators(h, h0, h1, h2, he, ht, h_G2);
    setupG_Generators(g, g0, g1, ge, g_G2);
    pairing(G, g, g_G2);

    // 사용자 및 키 초기화
    User user1, user2;
    user1.u.setByCSPRNG();
    user2.u.setByCSPRNG();

    std::string INFO = "send me a coin";

    Fp12::pow(user1.U, G, user1.u);
    Fp12::pow(user2.U, G, user2.u);

    alpha = 1;
    beta = 2;

    G2::mul(W, g_G2, alpha);
    G2::mul(X, h_G2, beta);

    // 측정 대상 user 준비
    User user;
    user.u.setByCSPRNG();
    Fp12::pow(user.U, G, user.u);

    if (account_est(user, g, g0, g1, g_G2, G)) {
        std::cout << "User 1 account created with signature A = " << user.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }

    // 시간 측정 시작
    std::ofstream csv("withdraw_coin_timing_origin.csv");
    csv << "Iteration,UserSide1_ns,BankSide_ns,UserSide2_ns\n";

    for (int i = 0; i < 10000; ++i) {
        long long user1_t = 0, bank_t = 0, user2_t = 0;

        bool success = withdraw_coin_time(user, h, h0, h1, h2, h_G2, G,
                                          user1_t, bank_t, user2_t);
        if (!success) {
            std::cerr << "Iteration " << i << ": withdraw_coin failed\n";
        }

        csv << i << "," << user1_t << "," << bank_t << "," << user2_t << "\n";
    }

    csv.close();
    sqlite3_close(db);
    std::cout << "[withdraw_coin] Measurement completed.\n";
}
void measure_payment_origin() {
    std::cout << "[payment] Measurement started...\n";
    if (!initDatabase(db)) return;
    initPairing(); // pairing 초기화

    // Generator 및 pairing 초기화
    G1 g, g0, g1, ge, h, h0, h1, h2, he, ht;
    G2 g_G2, h_G2;
    Fp12 G, H, H1;

    setupH_Generators(h, h0, h1, h2, he, ht, h_G2);
    setupG_Generators(g, g0, g1, ge, g_G2);
    pairing(G, g, g_G2);

    // 사용자 및 키 초기화
    User user1, user2;
    user1.u.setByCSPRNG();
    user2.u.setByCSPRNG();
    std::string INFO = "send me a coin";

    Fp12::pow(user1.U, G, user1.u);
    Fp12::pow(user2.U, G, user2.u);

    alpha = 1;
    beta  = 2;

    G2::mul(W, g_G2, alpha);
    G2::mul(X, h_G2, beta);

    // 계정 및 코인 발급
    if (account_est(user1, g, g0, g1, g_G2, G)) {
        std::cout << "User 1 account created with signature A = " << user1.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }

    if (account_est(user2, g, g0, g1, g_G2, G)) {
        std::cout << "User 2 account created with signature A = " << user2.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }

    if (withdraw_coin(user1, h, h0, h1, h2, h_G2, G)) {
        std::cout << "User 1 Coin created with signature B = " << user1.B.getStr(16) << std::endl;
    } else {
        std::cout << "Coin withdrawal failed." << std::endl;
    }

    // 시간 측정 시작
    std::ofstream csv("payment_timing_origin.csv");
    csv << "Iteration,Payee1_ns,Payer_ns,Payee2_ns\n";

    for (int i = 0; i < 10000; ++i) {
        long long payee1 = 0, payer_t = 0, payee2 = 0;

        bool success = Payment_time(user1, user2,
            g, g0, g1, g_G2, G,
            h, h0, h1, h2, h_G2,
            H, H1, INFO,
            payee1, payer_t, payee2);

        if (!success) {
            std::cerr << "Iteration " << i << ": Payment failed\n";
        }

        csv << i << "," << payee1 << "," << payer_t << "," << payee2 << "\n";
    }

    csv.close();
    sqlite3_close(db);
    std::cout << "[payment] Measurement completed.\n";
}
void measure_randomise_origin() {
    std::cout << "[randomise] Measurement started...\n";
    if (!initDatabase(db)) return;
    initPairing(); // pairing 초기화

    // Generator 및 pairing 초기화
    G1 g, g0, g1, ge, h, h0, h1, h2, he, ht;
    G2 g_G2, h_G2;
    Fp12 G, H, H1;

    setupH_Generators(h, h0, h1, h2, he, ht, h_G2);
    setupG_Generators(g, g0, g1, ge, g_G2);
    pairing(G, g, g_G2);

    // 사용자 및 키 초기화
    User user1, user2;
    user1.u.setByCSPRNG();
    user2.u.setByCSPRNG();
    std::string INFO = "send me a coin";

    Fp12::pow(user1.U, G, user1.u);
    Fp12::pow(user2.U, G, user2.u);

    alpha = 1;
    beta  = 2;
    G2::mul(W, g_G2, alpha);
    G2::mul(X, h_G2, beta);

    // 사전 단계: 계정 생성 및 코인 인출
    if (account_est(user1, g, g0, g1, g_G2, G)) {
        std::cout << "User 1 account created with signature A = " << user1.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }

    if (account_est(user2, g, g0, g1, g_G2, G)) {
        std::cout << "User 2 account created with signature A = " << user2.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }

    if (withdraw_coin(user1, h, h0, h1, h2, h_G2, G)) {
        std::cout << "User 1 Coin created with signature B = " << user1.B.getStr(16) << std::endl;
    } else {
        std::cout << "Coin withdrawal failed." << std::endl;
    }

    if (Payment(user1, user2, g, g0, g1, g_G2, G,
                h, h0, h1, h2, h_G2, H, H1, INFO)) {
        std::cout << "User 1 payment to User 2\n";
    } else {
        std::cout << "Payment failed.\n";
    }

    // 성능 측정
    std::ofstream csv("randomise_timing_origin.csv");
    csv << "Iteration,Payee1_ns,Bank1_ns,Payee2_ns,Bank2_ns,Payee3_ns\n";

    for (int i = 0; i < 10000; ++i) {
        long long t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0;
        bool success = randomise_time(user2,
                                      g, g0, g1, g_G2, G,
                                      h, h0, h1, h2, h_G2, H, H1,
                                      t1, t2, t3, t4, t5);
        if (!success) {
            std::cerr << "Iteration " << i << ": randomise failed\n";
        }

        csv << i << "," << t1 << "," << t2 << "," << t3 << "," << t4 << "," << t5 << "\n";
    }

    csv.close();
    sqlite3_close(db);
    std::cout << "[randomise] Measurement completed.\n";
}
void measure_finalise_origin() {
     std::cout << "[finalise] Measurement started...\n";
    if (!initDatabase(db)) return;
    initPairing(); // pairing 초기화

    // Generator 및 pairing 초기화
    G1 g, g0, g1, ge, h, h0, h1, h2, he, ht;
    G2 g_G2, h_G2;
    Fp12 G, H, H1;

    setupH_Generators(h, h0, h1, h2, he, ht, h_G2);
    setupG_Generators(g, g0, g1, ge, g_G2);
    pairing(G, g, g_G2);

    // 사용자 및 키 초기화
    User user1, user2;
    user1.u.setByCSPRNG();
    user2.u.setByCSPRNG();
    std::string INFO = "send me a coin";

    Fp12::pow(user1.U, G, user1.u);
    Fp12::pow(user2.U, G, user2.u);

    alpha = 1;
    beta = 2;
    G2::mul(W, g_G2, alpha);
    G2::mul(X, h_G2, beta);

    // 계정 생성
    if (account_est(user1, g, g0, g1, g_G2, G)) {
        std::cout << "User 1 account created with signature A = " << user1.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }
    if (account_est(user2, g, g0, g1, g_G2, G)) {
        std::cout << "User 2 account created with signature A = " << user2.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }

    // 코인 발행
    if (withdraw_coin(user1, h, h0, h1, h2, h_G2, G)) {
        std::cout << "User 1 Coin created with signature B = " << user1.B.getStr(16) << std::endl;
    } else {
        std::cout << "Coin withdrawal failed." << std::endl;
    }

    // 결제 수행
    if (Payment(user1, user2, g, g0, g1, g_G2, G,
                h, h0, h1, h2, h_G2, H, H1, INFO)) {
        std::cout << "User 1 payment to User 2\n";
    } else {
        std::cout << "Payment failed.\n";
    }

    // Finalise 성능 측정
    std::ofstream csv("finalise_timing.csv");
    csv << "Iteration,Payee1_ns,Bank1_ns\n";

    for (int i = 0; i < 10000; ++i) {
        long long t1 = 0, t2 = 0;
        bool success = finalise_time(user2,
                                     g, g0, g1, g_G2, G,
                                     h, h0, h1, h2, h_G2, H, H1,
                                     t1, t2);
        if (!success) {
            std::cerr << "Iteration " << i << ": finalise failed\n";
        }

        csv << i << "," << t1 << "," << t2 << "\n";
    }

    csv.close();
    sqlite3_close(db);
     std::cout << "[finalise] Measurement completed.\n";
}
int main() {

    std::cout << "===== Performance Measurement Started =====\n";

    measure_account_est_origin();
    measure_withdraw_coin_origin();
    measure_payment_origin();
    measure_randomise_origin();
    measure_finalise_origin();

    std::cout << "===== All Measurements Completed =====\n";
    return 0;
    /*
    if (!initDatabase(db)) return -1;
    initPairing(); // pairing 초기화

    G1 g,g0,g1,ge,h,h0,h1,h2,he,ht;
    G2 g_G2,h_G2;
    Fp12 G,H,H1;
    setupH_Generators(h, h0, h1, h2, he, ht, h_G2);
    setupG_Generators(g, g0, g1, ge, g_G2);
    pairing(G, g, g_G2);

    User user1,user2;
    user1.u.setByCSPRNG();
    user2.u.setByCSPRNG();
    std::string INFO="send me a coin";

    Fp12::pow(user1.U, G, user1.u);
    Fp12::pow(user2.U, G, user2.u);
    alpha=1;
    beta=2;

    G2::mul(W,g_G2,alpha);
    G2::mul(X,h_G2,beta);


    if (account_est(user1, g, g0, g1, g_G2, G)) {
        std::cout << "User 1 account created with signature A = " << user1.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }
    if (account_est(user2, g, g0, g1, g_G2, G)) {
        std::cout << "User 1 account created with signature A = " << user2.A.getStr(16) << std::endl;
    } else {
        std::cout << "Account establishment failed." << std::endl;
    }
    if (withdraw_coin(user1, h, h0, h1, h2, h_G2, G)) {
        std::cout << "User 1 Coin created with signature B = " << user1.B.getStr(16) << std::endl;
    } else {
        std::cout << "Coin withdrawal failed." << std::endl;
    }
    if (Payment(user1, user2, g, g0, g1, g_G2, G ,h, h0, h1, h2, h_G2, H, H1,INFO)) {
        std::cout << "User 1 payment to User 2\n";
    } else {
        std::cout << "Payment failed.\n";
    }
    if (finalise(user2, g, g0, g1, g_G2, G, h, h0, h1, h2 ,h_G2, H, H1)) {
        std::cout << "User 2 finalise coin\n";
    } else {
        std::cout << "Finalise failed.\n";
    }
    if (randomise( user2, g, g0, g1, g_G2, G ,h, h0, h1, h2, h_G2, H, H1)) {
        std::cout << "User 2 Randomise Coin\n";
    } else {
        std::cout << "Randomise failed.\n";
    }*/
    
    
    
    
    sqlite3_close(db);
    return 0;
}