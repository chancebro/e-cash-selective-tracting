#include <iostream>
#include "new.h"
#include <chrono>  // 이 헤더 추가
#include <fstream>

using namespace std::chrono;  // 네임스페이스 생략 위해
using namespace mcl::bn;

auto checkSetStr = [](auto& field, const char* valStr, const std::string& name, int base = 16) {
    try {
        using T = std::decay_t<decltype(field)>;
        if constexpr (std::is_same_v<T, mcl::bn::Fp12>) {
            std::stringstream ss(valStr);
            ss >> field;
            if (ss.fail()) throw cybozu::Exception("Fp12 load failed");
        } else {
            field.setStr(valStr, base);
        }
    } catch (const cybozu::Exception& e) {
        std::cerr << "❌ " << name << " setStr 실패\n값: " << valStr << "\n";
        std::cerr << "예외 메시지: " << e.what() << "\n";
        return false;
    }
    return true;
};


int main() {
    if (!initDatabase(db)) return -1;
    initPairing(); // pairing 초기화

    G1 g,g0,g1,gE,h,h0,h1,h2,hE,ht;
    G2 g_G2,h_G2;
    Fp12 G,GE,H,H1;
 
    setupH_Generators(h,  h0,  h1,  h2,hE,ht,h_G2);
    setupG_Generators(g,g0,g1,gE,g_G2);

    //개인키 값 초기화
    alpha=1;
    beta=2;
    x_GE=3;
    x_h=4;

    G2::mul(W,g_G2,alpha);
    G2::mul(X,h_G2,beta);

    pairing(G, g, g_G2);
    pairing(GE, g1, g_G2);
    pairing(H, h, h_G2);
    pairing(H1, h1, h_G2);
    
    User user1,user2,user3,user4,user5;

    user1.u.setByCSPRNG();
    user2.u.setByCSPRNG();
   
    std::string INFO;
    
    Fp12::pow(user1.U, G, user1.u);//generate User identifier U
    Fp12::pow(user2.U, G, user2.u);
    
    //insertUserU(db, 1, user1.U);
    //insertUserU(db, 2, user2.U);
    
    Fp12::pow(P,GE,x_GE);
    G1::mul(Q,hE,x_h);

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
    if (Payment(user1, user2, g, g0, g1, gE, g_G2, G,GE ,h, h0, h1, h2, hE, ht, h_G2, H, H1)) {
        std::cout << "User 1 payment to User 2\n";
    } else {
        std::cout << "Payment failed.\n";
    }
    std::ofstream csv("randomise_timing_proposal.csv");
    csv << "Iteration,Payee1_ns,Bank1_ns,Payee2_ns,Bank2_ns,Payee3_ns,Bank3_ns\n";

    for (int i = 0; i < 10000; ++i) {
        
        long long t1, t2, t3, t4, t5, t6;
        bool success = randomise_time(user2, g, g0, g1, g_G2, G, GE,
                                       h, h0, h1, h2, hE, ht, h_G2, H, H1,
                                       t1, t2, t3, t4, t5, t6);
        if (!success) std::cerr << "Randomise failed on iteration " << i << "\n";

        csv << i << "," << t1 << "," << t2 << "," << t3 << "," << t4 << "," << t5 << "," << t6 << "\n";
    }

    csv.close();
    sqlite3_close(db);
   /*//account gen of users
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
    
    if (Payment(user1, user2, g, g0, g1, gE, g_G2, G,GE ,h, h0, h1, h2, hE, ht, h_G2, H, H1)) {
        std::cout << "User 1 payment to User 2\n";
    } else {
        std::cout << "Payment failed.\n";
    }
     if (randomise(user2, g, g0, g1, g_G2, G, GE,h, h0, h1, h2,hE,ht ,h_G2, H, H1)) {
        std::cout << "User 2 randomise coin\n";
    } else {
        std::cout << "Randomise failed.\n";
    }
    if (finalise(user5, g, g0, g1, g_G2, G, GE,h, h0, h1, h2,hE,ht ,h_G2, H, H1)) {
        std::cout << "User 5 finalise coin\n";
    } else {
        std::cout << "Finalise failed.\n";
    }*/
    

   
}
