#include<bits/stdc++.h>
using namespace std;

// Tăng tốc độ I/O
void fast_io() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
}

int pop = 100;
const int T = 100;
int N;
int NE;
const double p = 1.0;
const double lenP = 5.0;
int Ne = 2;
// vector<vector<bool>> A; // Loại bỏ vì không cần thiết và tốn RAM
vector<vector<int>> x(pop + 1);
vector<vector<int>> e;
vector<int> xBest;
vector<int> d;
vector<vector<int>> dk(pop + 1);
vector<vector<int>> lk(pop + 1);
vector<int> pos(Ne);
random_device rd;
mt19937 gen(rd());
double rateMimicElist = 0.8;
int maLong = 0;
int macom = 0;

// Tính toán Delta Modularity nhanh: trả về lượng thay đổi Q nếu chuyển node u từ old_c sang new_c
double getDeltaQ(int u, int old_c, int new_c, const vector<int>& l_curr, const vector<int>& dk_curr, const vector<int>& lk_curr) {
    if (old_c == new_c) return 0.0;

    // Đếm số cạnh nối từ u tới cộng đồng cũ và mới
    // old_c != new_c (đã return ở trên) nên một label không thể vừa cũ vừa mới -> else if
    int k_u_old = 0;
    int k_u_new = 0;
    for (int v : e[u]) {
        int label_v = l_curr[v];
        if (label_v == old_c) k_u_old++;
        else if (label_v == new_c) k_u_new++;
    }

    double inv_m = 1.0 / (double)NE;
    double inv_m2 = 1.0 / (2.0 * NE);

    // Giá trị hiện tại của 2 cộng đồng
    double l_old = lk_curr[old_c];
    double d_old = dk_curr[old_c];
    double l_new = lk_curr[new_c];
    double d_new = dk_curr[new_c];

    // Đóng góp vào Q hiện tại của 2 cộng đồng này (dùng x*x thay cho pow(x,2.0))
    double a_old = d_old * inv_m2;
    double a_new = d_new * inv_m2;
    double q_old_term = (l_old * inv_m) - a_old * a_old;
    double q_new_term = (l_new * inv_m) - a_new * a_new;

    // Giá trị SAU khi chuyển
    // Rời old_c: giảm cạnh nội bộ k_u_old, giảm tổng bậc d[u]
    double l_old_next = l_old - k_u_old;
    double d_old_next = d_old - d[u];

    // Vào new_c: tăng cạnh nội bộ k_u_new, tăng tổng bậc d[u]
    double l_new_next = l_new + k_u_new;
    double d_new_next = d_new + d[u];

    double b_old = d_old_next * inv_m2;
    double b_new = d_new_next * inv_m2;
    double q_old_next_term = (l_old_next * inv_m) - b_old * b_old;
    double q_new_next_term = (l_new_next * inv_m) - b_new * b_new;

    // Delta = (Q_sau - Q_truoc)
    return (q_old_next_term + q_new_next_term) - (q_old_term + q_new_term);
}

// Hàm tính entropy và I được tối ưu hóa để tránh loop lồng quá lớn
double calI(const vector<int>& l1, const vector<int>& l2) {
    int s1 = 0, s2 = 0;
    for(int i=1; i<=N; i++) {
        if(l1[i] > s1) s1 = l1[i];
        if(l2[i] > s2) s2 = l2[i];
    }

    vector<int> c1_size(s1 + 1, 0);
    vector<int> c2_size(s2 + 1, 0);
    // Map: (nhãn l1 -> (Map: nhãn l2 -> số lượng))
    // Dùng vector phẳng để tăng tốc thay vì map
    vector<vector<int>> contingency(s1 + 1); 

    for (int i = 1; i <= N; i++) {
        c1_size[l1[i]]++;
        c2_size[l2[i]]++;
        
        // Chỉ lưu những node có label l1[i] tương ứng
        // Cách tối ưu: contingency[u] lưu danh sách các label l2 của các node thuộc l1=u
        contingency[l1[i]].push_back(l2[i]);
    }

    double I = 0;
    double N_double = (double)N;

    for (int i = 1; i <= s1; i++) {
        if (c1_size[i] == 0) continue;
        
        // Đếm tần suất các l2 trong nhóm l1=i
        // Để nhanh, ta sort hoặc dùng map nhỏ. Vì contingency[i] thường nhỏ so với N
        // Dùng map cục bộ hoặc sort
        sort(contingency[i].begin(), contingency[i].end());
        
        int current_l2 = -1;
        int count = 0;
        for(int val : contingency[i]){
            if(val != current_l2){
                if(count > 0){
                    double term = (double)count * N_double / ((double)c1_size[i] * (double)c2_size[current_l2]);
                    I += count * log(term);
                }
                current_l2 = val;
                count = 1;
            } else {
                count++;
            }
        }
        // Xử lý phần tử cuối
        if(count > 0){
             double term = (double)count * N_double / ((double)c1_size[i] * (double)c2_size[current_l2]);
             I += count * log(term);
        }
    }
    return I;
}

double calH(const vector<int>& l) {
    int s = 0;
    for(int i=1; i<=N; i++) if(l[i] > s) s = l[i];
    
    vector<int> c_size(s + 1, 0);
    for (int i = 1; i <= N; i++) c_size[l[i]]++;

    double H = 0;
    double N_double = (double)N;
    for (int i = 1; i <= s; i++) {
        if (c_size[i])
            H += (double)c_size[i] * log((double)c_size[i] / N_double);
    }
    return H;
}

double NMI(const vector<int>& l1, const vector<int>& l2) {
    return -2.0 * calI(l1, l2) / (calH(l1) + calH(l2));
}

double modularity(const vector<int>& dk, const vector<int>& lk) {
    double Q = 0;
    double inv_m = 1.0 / (double)NE;
    double inv_m2 = 1.0 / (2.0 * NE);
    int sz = lk.size();
    for (int i = 1; i < sz; i++) { // Chỉ duyệt các cộng đồng tồn tại
        if(dk[i] > 0 || lk[i] > 0) { // Check đơn giản
             double a = dk[i] * inv_m2;
             Q += (double)lk[i] * inv_m - a * a;
        }
    }
    return Q;
}

// Hàm thực hiện chuyển đổi (chỉ gọi khi đã confirm là tốt)
void applyTransfer(vector<int> &dk, vector<int> &lk, vector<int>& l, int u, int new_c) {
    int old_c = l[u];
    if(old_c == new_c) return;

    // Cập nhật dk
    dk[old_c] -= d[u];
    dk[new_c] += d[u];

    // Cập nhật lk
    // Phải duyệt qua neighbor để xem bao nhiêu cạnh nội bộ bị mất/thêm
    for (int v : e[u]) {
        if (l[v] == old_c) lk[old_c]--;
        if (l[v] == new_c) lk[new_c]++;
    }
    l[u] = new_c;
}

void LAR_rand(vector<vector<int>> &a) {
    for (int u = 1; u <= N; u++) {
        if (e[u].empty()) continue;
        uniform_int_distribution<int> disv(0, e[u].size() - 1);
        int v = e[u][disv(gen)];
        a[u].push_back(v);
        a[v].push_back(u);
    }
}

vector<int> decoding(const vector<vector<int>>& a) {
    vector<bool> dd(N + 1, false);
    vector<int> l(N + 1);
    int cnt = 0;

    for (int i = 1; i <= N; i++)
        if (!dd[i]) {
            ++cnt;
            queue<int> q;
            q.push(i);
            dd[i] = true; // Mark pushed
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                l[u] = cnt;
                for (int v : a[u])
                    if (!dd[v]) {
                        dd[v] = true;
                        q.push(v);
                    }
            }
        }
    return l;
}

void initialization() {
    for (int p = 1; p <= pop; p++) {
        vector<vector<int>> a(N + 1);
        LAR_rand(a);
        x[p] = decoding(a);

        int s = 0; 
        for(int i=1; i<=N; i++) if(x[p][i] > s) s = x[p][i];

        dk[p].assign(N + 1, 0); // Dùng assign nhanh hơn resize clear
        lk[p].assign(N + 1, 0);

        if((int)dk[p].size() <= s) dk[p].resize(s + 1, 0);
        if((int)lk[p].size() <= s) lk[p].resize(s + 1, 0);

        for (int u = 1; u <= N; u++) {
            int c = x[p][u];
            if(c >= dk[p].size()) { // Safety check
                dk[p].resize(c+1, 0); lk[p].resize(c+1, 0);
            }
            dk[p][c] += d[u];
            for (int v : e[u])
                if (x[p][u] == x[p][v] && u < v) {
                    lk[p][c]++;
                }
        }
    }
}

void movingToPrey(vector<int> &l, vector<int> &dk, vector<int> &lk, double k, const vector<int>& xTarget) {
    vector<int> pos_idx(N);
    iota(pos_idx.begin(), pos_idx.end(), 1);
    shuffle(pos_idx.begin(), pos_idx.end(), gen);

    int limit = min((int)k, N); 
    for (int i = 0; i < limit; i++) {
        int u = pos_idx[i];
        int target_label = xTarget[u];

        // --- SỬA LỖI SEGMENTATION FAULT TẠI ĐÂY ---
        // Nếu nhãn mục tiêu lớn hơn kích thước hiện có của dk/lk, hãy mở rộng chúng
        if (target_label >= (int)dk.size()) {
            dk.resize(target_label + 1, 0);
            lk.resize(target_label + 1, 0);
        }
        // ------------------------------------------

        applyTransfer(dk, lk, l, u, target_label);
    }
}

void randomWalk(vector<int> &l, vector<int> &dk, vector<int> &lk, double k) {
    vector<int> ranPop(pop);
    iota(ranPop.begin(), ranPop.end(), 1);
    shuffle(ranPop.begin(), ranPop.end(), gen);

    vector<int> randNode(N);
    iota(randNode.begin(), randNode.end(), 1);
    shuffle(randNode.begin(), randNode.end(), gen);

    int steps = (int)k;
    int idx = 0; // index for randNode
    
    // Logic vòng lặp i ngược từ lenP-1 về 0
    for (int i = lenP - 1; i >= 0; i--) {
        // BỔ SUNG AN TOÀN: Nếu population bị giảm xuống thấp hơn lenP (5), 
        // truy cập ranPop[i] sẽ lỗi nếu i >= pop.
        if (i >= (int)ranPop.size()) continue; 

        int kCommma = steps / (i + 1) + (int)(steps % (i + 1) > 0);
        while (kCommma-- && idx < N) {
            steps--;
            int u = randNode[idx++];
            
            // Lấy nhãn mục tiêu từ cá thể ngẫu nhiên P[i]
            int target_label = x[ranPop[i]][u];

            // --- SỬA LỖI SEGMENTATION FAULT TẠI ĐÂY ---
            // Kiểm tra và mở rộng bộ nhớ nếu nhãn mục tiêu vượt quá kích thước hiện tại
            if (target_label >= (int)dk.size()) {
                dk.resize(target_label + 1, 0);
                lk.resize(target_label + 1, 0);
            }
            // ------------------------------------------

            applyTransfer(dk, lk, l, u, target_label); 
        }
    }
}

void encirlingThePrey(vector<int>& l, vector<int> &dk, vector<int> &lk, double r, const vector<int>& xTarget) {
    int k = int(r * double(N));
    vector<int> pos_idx(N);
    iota(pos_idx.begin(), pos_idx.end(), 1);
    shuffle(pos_idx.begin(), pos_idx.end(), gen);

    int limit = min(k, N);
    for (int i = 0; i < limit; i++) {
        int u = pos_idx[i];
        int target_label = xTarget[u];

        // --- SỬA LỖI SEGMENTATION FAULT TẠI ĐÂY ---
        if (target_label >= (int)dk.size()) {
            dk.resize(target_label + 1, 0);
            lk.resize(target_label + 1, 0);
        }
        // ------------------------------------------

        applyTransfer(dk, lk, l, u, target_label);
    }
}

void mutation(vector<int> &l, vector<int> &dk, vector<int> &lk, double u_rate) {
    uniform_real_distribution<double> dis(0, 1);
    
    // Tìm max label hiện tại
    int S = 0;
    for(int val : l) if(val > S) S = val;

    const double inv_m  = 1.0 / (double)NE;
    const double inv_m2 = 1.0 / (2.0 * NE);
    // Đóng góp modularity của 1 cộng đồng (chỉ cần cho các cộng đồng bị ảnh hưởng)
    auto contrib = [&](int c) -> double {
        double a = dk[c] * inv_m2;
        return (double)lk[c] * inv_m - a * a;
    };

    // Buffer tái sử dụng để tránh cấp phát lại mỗi vòng lặp
    vector<int> nodes_to_move, old_labels, affected;

    for (int i = 1; i <= N; i++) {
        if (dis(gen) < u_rate) {
            // Thử tạo cộng đồng mới S+1
            int next_S = S + 1;
            // Lazy-resize: chỉ mở đủ tới next_S để các hàm modularity khác giữ vector nhỏ gọn
            if (next_S >= (int)dk.size()) {
                dk.resize(next_S + 1, 0);
                lk.resize(next_S + 1, 0);
            }

            double y = dis(gen);
            if (y < 0.5) {
                // Case 1: chỉ chuyển i sang cộng đồng mới — delta tính trực tiếp O(deg(i))
                double delta1 = getDeltaQ(i, l[i], next_S, l, dk, lk);
                if (delta1 > 0) {
                    applyTransfer(dk, lk, l, i, next_S);
                    S++;
                    macom = max(macom, S);
                }
            } else {
                // Case 2: chuyển i và toàn bộ hàng xóm sang cộng đồng mới.
                // Nhiều node chuyển cùng lúc, nhưng CHỈ các cộng đồng cũ của chúng và
                // cộng đồng mới next_S thay đổi -> tính delta Q trên đúng tập đó (O(deg))
                // thay vì quét toàn bộ vector như modularity() (O(N)). Kết quả tương đương.
                nodes_to_move.clear();
                nodes_to_move.push_back(i);
                for (int neighbor : e[i]) nodes_to_move.push_back(neighbor);

                // Tập cộng đồng bị ảnh hưởng = các nhãn cũ của node + next_S (loại trùng)
                affected.clear();
                for (int node : nodes_to_move) affected.push_back(l[node]);
                affected.push_back(next_S);
                sort(affected.begin(), affected.end());
                affected.erase(unique(affected.begin(), affected.end()), affected.end());

                double qBefore = 0;
                for (int c : affected) qBefore += contrib(c);

                old_labels.clear();
                for (int node : nodes_to_move) old_labels.push_back(l[node]);
                for (int node : nodes_to_move) applyTransfer(dk, lk, l, node, next_S);

                double qAfter = 0;
                for (int c : affected) qAfter += contrib(c);

                if (qAfter > qBefore) {
                    S++;
                    macom = max(macom, S);
                } else {
                    // Revert: l[node] hiện là next_S, chuyển trả về nhãn cũ
                    for (size_t k = 0; k < nodes_to_move.size(); k++) {
                        applyTransfer(dk, lk, l, nodes_to_move[k], old_labels[k]);
                    }
                }
            }
        }
    }
}

void caldklk(vector<int> &p,vector<int> &dk,vector<int> &lk){
    dk.assign(N+1,0); lk.assign(N+1,0);

    for (int u=1;u<=N;u++){
        dk[p[u]]+=d[u];
         
        for (int v:e[u])
            if (p[u]==p[v]&&u<v){
                ++lk[p[u]];
            }
    }
}

void consolidation(vector<int> &l, int l1, int l2){
    for (int i = 1; i <= N; i++)
        if (l[i] == l1)
            l[i] = l2;
}

void SecondaryCommunityConsolidation(vector<int> &l, vector<int> &dk, vector<int> &lk){
    // Thống kê cộng đồng + danh sách node theo nhãn (1 lần, O(N))
    unordered_map<int, int> cntNode;
    cntNode.reserve(N * 2);
    for (int i = 1; i <= N; i++) cntNode[l[i]]++;

    int numC = (int)cntNode.size();
    vector<pair<int,int>> decCommunities; // (nhãn, kích thước) — giảm dần theo size
    decCommunities.reserve(numC);
    for (auto &kv : cntNode) decCommunities.push_back({kv.first, kv.second});
    sort(decCommunities.begin(), decCommunities.end(),
         [](const pair<int,int>& a, const pair<int,int>& b){ return a.second > b.second; });

    // Danh sách node của từng cộng đồng (được cập nhật khi gộp)
    unordered_map<int, vector<int>> nodesInComm;
    nodesInComm.reserve(numC * 2);
    for (int i = 1; i <= N; i++) nodesInComm[l[i]].push_back(i);

    const double inv_m  = 1.0 / (double)NE;
    const double inv_m2 = 1.0 / (2.0 * NE);
    auto contrib = [&](int c) -> double {
        double a = dk[c] * inv_m2;
        return (double)lk[c] * inv_m - a * a;
    };

    int i = numC - 1;
    while (i > 0){
        int A = decCommunities[i].first;
        int j = 0;
        bool check = false;
        while (i > j){
            int B = decCommunities[j].first;

            // Số cạnh nối giữa cộng đồng A và B (duyệt node của A — thường là cộng đồng nhỏ)
            int cross = 0;
            for (int u : nodesInComm[A])
                for (int v : e[u])
                    if (l[v] == B) cross++;

            // Gộp A vào B chỉ thay đổi đóng góp Q của A và B -> delta tính trực tiếp
            double before = contrib(A) + contrib(B);
            int newLkB = lk[A] + lk[B] + cross;
            int newDkB = dk[A] + dk[B];
            double aB = newDkB * inv_m2;
            double after = (double)newLkB * inv_m - aB * aB; // A trở thành rỗng -> đóng góp 0

            if (after > before){
                // Thực hiện gộp A -> B
                for (int u : nodesInComm[A]){
                    l[u] = B;
                    nodesInComm[B].push_back(u);
                }
                nodesInComm[A].clear();
                lk[B] = newLkB; dk[B] = newDkB;
                lk[A] = 0;      dk[A] = 0;
                --i;
                check = true;
                break;
            }
            ++j;
        }
        if (!check) --i;
    }
}

void boudaryNodeAdjustment(vector<int> &l, vector<int> &dk, vector<int> &lk, double rate) {
    uniform_real_distribution<double> dis(0.0, 1.0);
    
    // Tối ưu: Không cần copy vector tmpl, dktmp, lktmp liên tục
    for (int i = 1; i <= N; i++) {
        if (dis(gen) > rate) continue;
        
        // Chỉ xét những neighbor khác cộng đồng
        int best_c = l[i];
        double best_delta = 0;
        bool found = false;

        // Tìm cộng đồng hàng xóm tốt nhất để chuyển sang (Greedy optimization)
        // Code cũ: chuyển sang neighbor ĐẦU TIÊN thấy khác và check modularity
        // Để giữ logic cũ:
        for (int neighbor : e[i]) {
            if (l[i] != l[neighbor]) {
                 // Tính thử Delta nếu chuyển i sang l[neighbor]
                 double delta = getDeltaQ(i, l[i], l[neighbor], l, dk, lk);
                 if(delta > 0) {
                     // Nếu tốt hơn thì chuyển luôn (Greedy) giống logic cũ (nhưng nhanh hơn)
                     applyTransfer(dk, lk, l, i, l[neighbor]);
                     // Break để sang node tiếp theo (giống logic code cũ là thử chuyển và nếu tốt giữ lại)
                     // Code gốc: lặp neighbor, nếu khác -> chuyển thử -> nếu tệ -> revert.
                     // Code cũ KHÔNG break, nó thử hết các neighbor? 
                     // Code cũ: `if (l[i]!=l[neighbor] && !dd[l[neighbor]])` -> dd chỉ đánh dấu label đã thử
                     // Nên logic là thử tất cả các label hàng xóm DUY NHẤT.
                     
                     // Ở đây ta đơn giản hóa: Nếu chuyển tốt thì giữ, và tiếp tục vòng lặp neighbor
                     // (tức là có thể chuyển nhiều lần cho 1 node trong 1 lần gọi hàm, hoặc chuyển sang A rồi sang B)
                     // Tuy nhiên, việc tính Delta dựa trên state hiện tại là chính xác.
                 }
            }
        }
    }
}

void EPD() {
    if (x.size() <= 10) return;

    vector<pair<double, int>> modularityValues;
    modularityValues.reserve(pop);
    for (int i = 1; i <= pop; i++) {
        double modValue = modularity(dk[i], lk[i]);
        modularityValues.push_back({modValue, i});
    }

    sort(modularityValues.rbegin(), modularityValues.rend());

    // Sắp xếp lại x, dk, lk theo thứ tự tốt nhất
    // Tối ưu: Dùng vector hoán vị thay vì copy data mảng lớn nếu có thể, 
    // nhưng ở đây copy vector<int> không quá chậm với pop=100.
    vector<vector<int>> newX(pop + 1);
    vector<vector<int>> newDk(pop + 1);
    vector<vector<int>> newLk(pop + 1);

    for (int i = 0; i < pop; i++) {
        int idx = modularityValues[i].second;
        newX[i + 1] = move(x[idx]); // Dùng move để tránh deep copy
        newDk[i + 1] = move(dk[idx]);
        newLk[i + 1] = move(lk[idx]);
    }
    x = move(newX);
    dk = move(newDk);
    lk = move(newLk);

    // Find EL
    // Chỉ giữ logic tìm kiếm, tối ưu loop
    // Logic của bạn là loại bỏ các cá thể giống nhau (NMI=1) hoặc kém.
    // Logic này giữ nguyên.
    
    // Lưu ý: Sau khi sort, x[1]..x[Ne] là tốt nhất.
    for(int i=0; i<Ne; i++) {
        if (i + 1 <= pop) pos[i] = i + 1;
        else pos[i] = 1; // Fallback về 1 nếu pop quá nhỏ
    }
    // ----------------------------------------------------

    int cnt = 0;
    // Đoạn logic tìm kiếm Elite giữ nguyên, nhưng thay đổi cách gán pos
    // Lưu ý: pos đã có giá trị mặc định, ta chỉ update nếu tìm thấy cái tốt hơn/khác biệt hơn
    
    // Logic cũ của bạn hơi phức tạp và dễ sai index khi erase. 
    // Hãy dùng logic đơn giản hóa này để tìm Elite đa dạng:
    
    vector<int> current_elites;
    if (pop >= 1) current_elites.push_back(1); // Luôn lấy con tốt nhất
    
    for (int i = 2; i <= pop && current_elites.size() < Ne; i++) {
        bool distinct = true;
        for (int elite_idx : current_elites) {
            if (NMI(x[i], x[elite_idx]) > 0.99) { // Nếu quá giống elite đã chọn
                distinct = false;
                break;
            }
        }
        if (distinct) {
            current_elites.push_back(i);
        }
    }
    
    // Cập nhật lại pos từ danh sách đã tìm được
    for(int i=0; i<Ne; i++) {
        if (i < current_elites.size()) pos[i] = current_elites[i];
        else pos[i] = 1; // Fallback an toàn
    }

    // --- Phần xóa phần tử trùng lặp (Real EPD removal) ---
    // Cần cẩn thận: khi xóa, index của pos sẽ bị sai lệch.
    // Để an toàn và đơn giản, ta KHÔNG xóa x trong hàm EPD nữa (hoặc chỉ xóa đuôi).
    // Việc xóa phần tử giữa vector làm thay đổi index của pos[0], pos[1] rất nguy hiểm.
    
    // NẾU BẠN MUỐN GIỮ LOGIC XÓA:
    // Hãy xóa từ dưới lên (từ pop về 1) để hạn chế ảnh hưởng index nhỏ.
    
    double N_nor = pop - (pop / 2 + 1) + 1;
    uniform_real_distribution<double> dis(0, 1);
    
    // Duyệt ngược để erase an toàn hơn
    for (int i = pop; i >= pop / 2 + 1; i--) {
        bool isElite = false;
        for(int k=0; k<Ne; k++) if(i == pos[k]) isElite = true;
        
        if (!isElite) {
            double C = 1.0 - exp(-double(i) / N_nor);
            if (dis(gen) <= C) {
                // Không cần chỉnh sửa pos vì ta đang xóa ở nửa dưới (index lớn),
                // còn pos thường nằm ở nửa trên (index nhỏ - fitness cao).
                // Tuy nhiên vẫn cần check:
                bool affect_pos = false;
                for(int k=0; k<Ne; k++) if(pos[k] > i) affect_pos = true;
                
                if (!affect_pos) { // Chỉ xóa nếu không ảnh hưởng vị trí elite
                    x.erase(x.begin() + i);
                    dk.erase(dk.begin() + i);
                    lk.erase(lk.begin() + i);
                    pop--;
                }
            }
        }
    }
    
    if (Ne > 0 && pos[Ne-1] > 0 && pos[Ne-1] <= pop)
        maLong = NMI(x[pos[Ne - 1]], xBest);
}

void updateLocation(vector<int> &l, int t, vector<int> &dk, vector<int> &lk) {
    uniform_real_distribution<double> dis(0, 1);
    double alpha = dis(gen), beta = dis(gen);
    
    const vector<int>* target = &xBest; // Mặc định an toàn nhất là xBest

    // Helper lambda để lấy elite an toàn
    auto getSafeElite = [&](int idx_in_pos) -> const vector<int>* {
        int idx = pos[idx_in_pos];
        if (idx >= 1 && idx <= pop) return &x[idx];
        return &xBest; // Fallback nếu pos lỗi
    };

    if (alpha < 0.5) {
        double k = p * double(t) * double(N) / double(T);
        if (beta < 0.5) {
            if (dis(gen) < rateMimicElist) target = &xBest;
            else {
                if (dis(gen) < 0.5) target = getSafeElite(0);
                else target = getSafeElite(1);
            }
            movingToPrey(l, dk, lk, k, *target);
        } else {
            randomWalk(l, dk, lk, k);
        }
    } else {
        uniform_real_distribution<double> disl(-1, 1);
        double ll = disl(gen);
        double r = abs(cos(2 * M_PI * ll));

        if (dis(gen) < rateMimicElist) target = &xBest;
        else {
            if (dis(gen) < 0.5) target = getSafeElite(0);
            else target = getSafeElite(1);
        }
        encirlingThePrey(l, dk, lk, r, *target);
    }
}

void standardization(vector<int> &p) {
    map<int, int> mp;
    int cnt = 0;
    for (int i = 1; i <= N; i++)
        if (!mp[p[i]])
            mp[p[i]] = ++cnt;
    for (int i = 1; i <= N; i++)
        p[i] = mp[p[i]];
}

void EP_WOCD() {
    initialization();
    double ans = -1e9;
    
    // Tìm best ban đầu
    for (int i = 1; i <= pop; i++) {
        double val = modularity(dk[i], lk[i]);
        if (val > ans) {
            ans = val;
            xBest = x[i];
        }
    }

    for (int t = 1; t <= T; t++) {
        double ib = -1e9;
        
        // Loop chính
        for (int p_idx = 1; p_idx <= pop; p_idx++) {
            double rateLS = 1.0; 
            double rateMu = 0.3;
            bool check = (p_idx > Ne);

            if (check) updateLocation(x[p_idx], t, dk[p_idx], lk[p_idx]);
            
            mutation(x[p_idx], dk[p_idx], lk[p_idx], rateMu);
            boudaryNodeAdjustment(x[p_idx], dk[p_idx], lk[p_idx], rateLS);
        }

        // Cập nhật Best Global sau mỗi vòng lặp
        for (int i = 1; i <= pop; i++) {
            double val = modularity(dk[i], lk[i]);
            if (val > ans) {
                ans = val;
                xBest = x[i];
            }
            if (val > ib) ib = val;
        }

        EPD();
    }
    
    vector<int> dkBest, lkBest;
    dkBest.assign(N + 1, 0);
    lkBest.assign(N + 1, 0);
    for (int u = 1; u <= N; u++) {
        int c = xBest[u];
        if(c >= dkBest.size()) { // Safety check
            dkBest.resize(c+1, 0); lkBest.resize(c+1, 0);
        }
        dkBest[c] += d[u];
        for (int v : e[u])
            if (xBest[u] == xBest[v] && u < v) {
                lkBest[c]++;
            }
    }
    SecondaryCommunityConsolidation(xBest, dkBest, lkBest);
    ans=modularity(dkBest,lkBest);

    for (int i=1;i<=pop;i++){
        SecondaryCommunityConsolidation(x[i],dk[i],lk[i]);
        if (modularity(dk[i],lk[i])>ans){
            ans=modularity(dk[i],lk[i]);
            xBest=x[i];
        }
    }
    // Modularity ra stderr để stdout chỉ chứa nhãn cộng đồng (cho thư viện Python đọc)
    cerr << ans << "\n";
}

// Dùng làm thư viện: đọc đồ thị từ stdin, in nhãn cộng đồng của từng node ra stdout.
//
// Định dạng đầu vào (stdin):
//   N NE
//   u v        (NE dòng, node đánh số 1..N, vô hướng)
// Tham số dòng lệnh:
//   argv[1] (tuỳ chọn) = seed cho bộ sinh ngẫu nhiên (để tái lập kết quả)
// Đầu ra (stdout):
//   N dòng, dòng thứ i là nhãn cộng đồng (đã chuẩn hoá 1..k) của node i.
int main(int argc, char** argv) {
    fast_io();

    // Seed: nếu được truyền vào thì tái lập được; nếu không thì ngẫu nhiên.
    if (argc > 1) {
        gen.seed((unsigned)strtoull(argv[1], nullptr, 10));
    }

    if (!(cin >> N >> NE)) return 0;

    d.assign(N + 1, 0);
    e.assign(N + 1, {});

    for (int i = 1; i <= NE; i++) {
        int u, v;
        cin >> u >> v;
        e[u].push_back(v);
        e[v].push_back(u);
        d[u]++;
        d[v]++;
    }

    EP_WOCD();

    // Chuẩn hoá nhãn về 1..k rồi xuất nhãn của từng node.
    standardization(xBest);
    string out;
    out.reserve(N * 4);
    for (int i = 1; i <= N; i++) {
        out += to_string(xBest[i]);
        out += '\n';
    }
    cout << out;
    return 0;
}