#include <iostream>
#include <sstream>
#include <sqlite3.h>
#include "mcl/bn256.hpp"
#include "new.h"

using namespace mcl::bn;


int forwardTracing(sqlite3* db, int start_ts_num) {
    Fp12 S, userID_C1_payer, userID_C2_payer, userID_C1_payee, userID_C2_payee;
    G1 T, backward_C1, backward_C2, forward_C1, forward_C2, bank_B;

    while (start_ts_num != -1) {
        const char* sql_template = "SELECT S, T, backward_C1, backward_C2, bank_B, forward_C1, forward_C2, "
                                   "userID_payer_C1, userID_payer_C2, userID_payee_C1, userID_payee_C2 "
                                   "FROM spk_bundle WHERE Ts_num = ?;";
        sqlite3_stmt* stmt;
        if (sqlite3_prepare_v2(db, sql_template, -1, &stmt, nullptr) != SQLITE_OK) {
            std::cerr << "❌ SELECT 구문 준비 실패\n";
            return -1;
        }
        sqlite3_bind_int(stmt, 1, start_ts_num);

        bool finalised = false;
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            const char* bank_b_raw = (const char*)sqlite3_column_text(stmt, 4);
            const char* fwd1_raw = (const char*)sqlite3_column_text(stmt, 5);
            const char* fwd2_raw = (const char*)sqlite3_column_text(stmt, 6);

            if (!bank_b_raw || !fwd1_raw || !fwd2_raw || 
                std::string(bank_b_raw) == "none" || 
                std::string(fwd1_raw) == "none" || 
                std::string(fwd2_raw) == "none") {
                finalised = true;
            }

            try {
                std::stringstream ssS((const char*)sqlite3_column_text(stmt, 0)); ssS >> S;
                T.setStr((const char*)sqlite3_column_text(stmt, 1), 16);
                backward_C1.setStr((const char*)sqlite3_column_text(stmt, 2), 16);
                backward_C2.setStr((const char*)sqlite3_column_text(stmt, 3), 16);
                if (!finalised) {
                    bank_B.setStr(bank_b_raw, 16);
                    forward_C1.setStr(fwd1_raw, 16);
                    forward_C2.setStr(fwd2_raw, 16);
                }
                std::stringstream((const char*)sqlite3_column_text(stmt, 7)) >> userID_C1_payer;
                std::stringstream((const char*)sqlite3_column_text(stmt, 8)) >> userID_C2_payer;
                std::stringstream((const char*)sqlite3_column_text(stmt, 9)) >> userID_C1_payee;
                std::stringstream((const char*)sqlite3_column_text(stmt, 10)) >> userID_C2_payee;
            } catch (const cybozu::Exception& e) {
                std::cerr << "❌ 필드 파싱 실패: " << e.what() << "\n";
                sqlite3_finalize(stmt);
                return -1;
            }
        } else {
            std::cout << "❌ 해당 Ts_num을 찾을 수 없습니다.\n";
            sqlite3_finalize(stmt);
            return -1;
        }
        sqlite3_finalize(stmt);

        // 복호화: recovered_U = C2 / C1^x_GE
        Fp12 temp_pow, recovered_U_payer;
        Fp12::pow(temp_pow, userID_C1_payer, x_GE);
        Fp12::inv(temp_pow, temp_pow);
        Fp12::mul(recovered_U_payer, userID_C2_payer, temp_pow);

        Fp12 recovered_U_payee;
        Fp12::pow(temp_pow, userID_C1_payee, x_GE);
        Fp12::inv(temp_pow, temp_pow);
        Fp12::mul(recovered_U_payee, userID_C2_payee, temp_pow);

        std::string recoveredStr_payer = recovered_U_payer.getStr(16);
        std::string recoveredStr_payee = recovered_U_payee.getStr(16);

        int payer_id = -1, payee_id = -1;
        const char* sql_U_table = "SELECT id, U FROM users;";
        sqlite3_stmt* stmt_u;
        if (sqlite3_prepare_v2(db, sql_U_table, -1, &stmt_u, nullptr) != SQLITE_OK) {
            std::cerr << "❌ U 테이블 SELECT 준비 실패\n";
            return -1;
        }

        while (sqlite3_step(stmt_u) == SQLITE_ROW) {
            int row_id = sqlite3_column_int(stmt_u, 0);
            std::string uStr = reinterpret_cast<const char*>(sqlite3_column_text(stmt_u, 1));
            if (uStr == recoveredStr_payer) payer_id = row_id;
            if (uStr == recoveredStr_payee) payee_id = row_id;
        }
        sqlite3_finalize(stmt_u);

        std::cout << "👤 Payer ID: " << (payer_id != -1 ? std::to_string(payer_id) : "❌ 매칭 실패") << "\n";
        std::cout << "👤 Payee ID: " << (payee_id != -1 ? std::to_string(payee_id) : "❌ 매칭 실패") << "\n";

        if (finalised) {
            std::cout << "✅ Finalised 상태 감지됨 - Tracing 종료\n";
            return 0;
        }

        // Step 1: T = forward_C2 - x_h * forward_C1 + bank_B
        G1 temp_check_T, fC1_xh;
        G1::mul(fC1_xh, forward_C1, x_h);
        G1::sub(temp_check_T, forward_C2, fC1_xh);
        G1::add(temp_check_T, temp_check_T, bank_B);

        const char* sql_next = "SELECT Ts_num FROM spk_bundle WHERE T = ?;";//이렇게 T를 가진 스랜잭션을 찾음음
        sqlite3_stmt* stmt_next;
        if (sqlite3_prepare_v2(db, sql_next, -1, &stmt_next, nullptr) != SQLITE_OK) {
            std::cerr << "❌ 다음 T SELECT 준비 실패\n";
            return -1;
        }

        std::string tStr = temp_check_T.getStr(16);
        sqlite3_bind_text(stmt_next, 1, tStr.c_str(), -1, SQLITE_STATIC);

        start_ts_num = -1;
        if (sqlite3_step(stmt_next) == SQLITE_ROW) {
            start_ts_num = sqlite3_column_int(stmt_next, 0);
            std::cout << "🔁 다음 연결된 Ts_num: " << start_ts_num << "\n";
        } else {
            std::cout << "✅ Tracing 종료: 더 이상 연결된 T 없음\n";
        }
        sqlite3_finalize(stmt_next);
    }
    return 0;
}

int BackwardTracing(sqlite3* db, int start_ts_num) {
    Fp12 S, userID_C1_payer, userID_C2_payer, userID_C1_payee, userID_C2_payee;
    G1 T, backward_C1, backward_C2, forward_C1, forward_C2, bank_B;

    while (start_ts_num != -1) {
        const char* sql_template = "SELECT S, T, backward_C1, backward_C2, bank_B, forward_C1, forward_C2, "
                                   "userID_payer_C1, userID_payer_C2, userID_payee_C1, userID_payee_C2 "
                                   "FROM spk_bundle WHERE Ts_num = ?;";
        sqlite3_stmt* stmt;
        if (sqlite3_prepare_v2(db, sql_template, -1, &stmt, nullptr) != SQLITE_OK) {
            std::cerr << "❌ SELECT 구문 준비 실패\n";
            return -1;
        }
        sqlite3_bind_int(stmt, 1, start_ts_num);

        if (sqlite3_step(stmt) == SQLITE_ROW) {
            try {
                std::stringstream ssS((const char*)sqlite3_column_text(stmt, 0)); ssS >> S;
                T.setStr((const char*)sqlite3_column_text(stmt, 1), 16);
                backward_C1.setStr((const char*)sqlite3_column_text(stmt, 2), 16);
                backward_C2.setStr((const char*)sqlite3_column_text(stmt, 3), 16);

                const char* bank_b_raw = (const char*)sqlite3_column_text(stmt, 4);
                const char* fwd1_raw = (const char*)sqlite3_column_text(stmt, 5);
                const char* fwd2_raw = (const char*)sqlite3_column_text(stmt, 6);

                if (bank_b_raw && std::string(bank_b_raw) != "none") {
                    bank_B.setStr(bank_b_raw, 16);
                }
                if (fwd1_raw && std::string(fwd1_raw) != "none") {
                    forward_C1.setStr(fwd1_raw, 16);
                }
                if (fwd2_raw && std::string(fwd2_raw) != "none") {
                    forward_C2.setStr(fwd2_raw, 16);
                }

                std::stringstream((const char*)sqlite3_column_text(stmt, 7)) >> userID_C1_payer;
                std::stringstream((const char*)sqlite3_column_text(stmt, 8)) >> userID_C2_payer;
                std::stringstream((const char*)sqlite3_column_text(stmt, 9)) >> userID_C1_payee;
                std::stringstream((const char*)sqlite3_column_text(stmt, 10)) >> userID_C2_payee;
            } catch (const cybozu::Exception& e) {
                std::cerr << "❌ 필드 파싱 실패: " << e.what() << "\n";
                sqlite3_finalize(stmt);
                return -1;
            }
        } else {
            std::cout << "❌ 해당 Ts_num을 찾을 수 없습니다.\n";
            sqlite3_finalize(stmt);
            return -1;
        }
        sqlite3_finalize(stmt);

        //std::cout << "\n🔁 Backward Tracing 시작 Ts_num: " << start_ts_num << "\n";

        // Payer / Payee ID 복호화
        Fp12 temp_pow_fp12, recovered_U_payer;
        Fp12::pow(temp_pow_fp12, userID_C1_payer, x_GE);
        Fp12::inv(temp_pow_fp12, temp_pow_fp12);
        Fp12::mul(recovered_U_payer, userID_C2_payer, temp_pow_fp12);

        Fp12 recovered_U_payee;
        Fp12::pow(temp_pow_fp12, userID_C1_payee, x_GE);
        Fp12::inv(temp_pow_fp12, temp_pow_fp12);
        Fp12::mul(recovered_U_payee, userID_C2_payee, temp_pow_fp12);

        std::string recoveredStr_payer = recovered_U_payer.getStr(16);
        std::string recoveredStr_payee = recovered_U_payee.getStr(16);

        int payer_id = -1, payee_id = -1;
        const char* sql_U_table = "SELECT id, U FROM users;";
        sqlite3_stmt* stmt_u;
        if (sqlite3_prepare_v2(db, sql_U_table, -1, &stmt_u, nullptr) != SQLITE_OK) {
            std::cerr << "❌ U 테이블 SELECT 준비 실패\n";
            return -1;
        }

        while (sqlite3_step(stmt_u) == SQLITE_ROW) {
            int row_id = sqlite3_column_int(stmt_u, 0);
            std::string uStr = reinterpret_cast<const char*>(sqlite3_column_text(stmt_u, 1));
            if (uStr == recoveredStr_payer) payer_id = row_id;
            if (uStr == recoveredStr_payee) payee_id = row_id;
        }
        sqlite3_finalize(stmt_u);

        std::cout << "👤 Payer ID: " << (payer_id != -1 ? std::to_string(payer_id) : "❌ 매칭 실패") << "\n";
        std::cout << "👤 Payee ID: " << (payee_id != -1 ? std::to_string(payee_id) : "❌ 매칭 실패") << "\n";

        // 복호화: recovered_ht_t = backward_C2 - x_h * backward_C1
        G1 temp_pow, neg_pow, recovered_ht_t;
        G1::mul(temp_pow, backward_C1, x_h);
        G1::neg(neg_pow, temp_pow);
        G1::add(recovered_ht_t, backward_C2, neg_pow);

        // target_bank_B = T - ht^t
        G1 inv_ht_t;
        G1::neg(inv_ht_t, recovered_ht_t);
        G1 target_bank_B;
        G1::add(target_bank_B, T, inv_ht_t);

        // DB에서 target_bank_B를 가진 행 찾기
        const char* sql_find = "SELECT Ts_num FROM spk_bundle WHERE bank_B = ?;";
        sqlite3_stmt* stmt_next;
        if (sqlite3_prepare_v2(db, sql_find, -1, &stmt_next, nullptr) != SQLITE_OK) {
            std::cerr << "❌ 다음 bank_B SELECT 준비 실패\n";
            return -1;
        }
        std::string bStr = target_bank_B.getStr(16);
        sqlite3_bind_text(stmt_next, 1, bStr.c_str(), -1, SQLITE_STATIC);

        int next_ts = -1;
        if (sqlite3_step(stmt_next) == SQLITE_ROW) {
            next_ts = sqlite3_column_int(stmt_next, 0);
        }
        sqlite3_finalize(stmt_next);

        if (next_ts != -1) {
            std::cout << "🔁 다음 연결된 Ts_num: " << next_ts << "\n";
            start_ts_num = next_ts;
        } else {
            std::cout << "✅ Tracing 종료: 더 이상 연결된 bank_B 없음\n";
            break;
        }
    }
    return 0;
}

int main() {
    initPairing();
    G1 g,g0,g1,gE,h,h0,h1,h2,hE,ht;
    G2 g_G2,h_G2;
    Fp12 G,GE,H,H1;

    setupH_Generators(h, h0, h1, h2, hE, ht, h_G2);
    setupG_Generators(g, g0, g1, gE, g_G2);

    alpha = 1; beta = 2; x_GE = 3; x_h = 4;
    G2::mul(W, g_G2, alpha);
    G2::mul(X, h_G2, beta);

    pairing(G, g, g_G2);
    pairing(GE, g1, g_G2);
    pairing(H, h, h_G2);
    pairing(H1, h1, h_G2);

    Fp12::pow(P, GE, x_GE);
    G1::mul(Q, hE, x_h);

    sqlite3* db;
    if (sqlite3_open("Tracing_sample.db", &db) != SQLITE_OK) {
        std::cerr << "❌ DB 열기 실패\n";
        return 1;
    }

    int startTs_forward = 1; // 최초 감사 대상 트랜잭션 Ts_num 입력
    std::cout << "🔁 Forward Tracing 시작 Ts_num: " << startTs_forward<< "\n";
    forwardTracing(db, startTs_forward);

    int startTs_backward = 4; // 최초 감사 대상 트랜잭션 Ts_num 입력
    std::cout << "🔁 Backward Tracing 시작 Ts_num: " << startTs_backward<< "\n";
    BackwardTracing(db, startTs_backward);
    sqlite3_close(db);
    return 0;
}
