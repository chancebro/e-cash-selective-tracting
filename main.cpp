#include <iostream>
#include <sqlite3.h>
#include <mcl/bn.hpp>
#include "tracing.h"
#include <chrono>  // 이 헤더 추가
#include <fstream>

using namespace std::chrono;  // 네임스페이스 생략 위해
using namespace mcl::bn;



int main() {
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
    std::string INFO;

    Fp12::pow(user1.U, G, user1.u);
    Fp12::pow(user2.U, G, user2.u);

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
    
    std::ofstream csv("randomise_timing.csv");
    csv << "Iteration,Payee1_ns,Bank1_ns,Payee2_ns,Bank2_ns,Payee3_ns\n";
   
    for (int i = 0; i < 10000; ++i) {

        long long t1, t2, t3, t4, t5;
        bool success = randomise_time(user2, g, g0, g1, g_G2, G,
                                       h, h0, h1, h2, h_G2, H, H1,
                                       t1, t2, t3, t4, t5);
        if (!success) {
            std::cerr << "Iteration " << i << ": randomise failed\n";
        }

        csv << i << "," << t1 << "," << t2 << "," << t3 << "," << t4 << "," << t5 << "\n";
    }

    csv.close();
    
    sqlite3_close(db);
    return 0;
}