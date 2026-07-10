clear
rng(0)

Eb_No_db =15.0;

T=23; %max decoder iteration

eta=12; 
w=30;

flip_num=2;

p=3; % number of bits per symbol
q = 2^p; 

load('arith_8.mat');
load('nbldpc_45_2_3_gf8.mat'); 

[G, K, pivotCols, freeCols, H_rref] = gf8_generator_from_H(h, add_mat, mul_mat, div_mat);


%h = full(h); % H matrix
N = size(h,2); %length of codeword
M = size(h,1); %number of parity checks
%K=N-M; % msg length
R=K/N; % code rate

[M, ~] = size(h);
str_cn_vn = cell(M,1);

for i = 1:M
    str_cn_vn{i} = find(h(i,:) ~= 0);
end

CN_lst = str_cn_vn;

info_seq=randi([0 q-1], 1, K);

code_seq = gf_vec_mat_mul(info_seq, G, add_mat, mul_mat);

Syndromes = decod_prod(code_seq, h, CN_lst, mul_mat, add_mat);

%QAM-8 mapping
[qam8, qam_binary_map] = generate_qam8_map();

avg_pow = qam8' * qam8 / 8; %sum of squared magnitudes of all symbols(the total power)/q.
nrm_fct=sqrt(avg_pow); %normalization factor
gf8 = (0:q-1); %GF field symbols
%alph_bin =  logical(fliplr(dec2bin(gf8, p) - 48)); % symbols in binary

%initializing vector size

y = zeros(1, N);
hard_d_cmplx = zeros(1, N);
hard_d_gf8  = zeros(1, N);


FE=zeros(length(Eb_No_db),1);
genFrame=zeros(length(Eb_No_db),1);
iters_cnr=zeros(length(Eb_No_db),1);
BE_Coded=zeros(length(Eb_No_db),1);
BE_unCoded=zeros(length(Eb_No_db),1);


targetFE=150; %maximum FER to be observed
max_gen=1e6; % maximum number of frame to be generated


fail_log = struct([]);
fail_count = 0;

for i = 1:length(Eb_No_db)

     while(FE(i) < targetFE && genFrame(i)<max_gen)

        genFrame(i)=genFrame(i)+1;

        c(1,1:N) = qam8(code_seq'+1,1); % codeword in complex
        avg_symbol_energy = 1;

        Eb_No_linear(i)= 10.^(Eb_No_db(i)/10);
        No = avg_symbol_energy ./ (p * Eb_No_linear(i)*R); %noise spectral density
        sigma0 = sqrt(No/2)*nrm_fct ; %noise standard deviation
        nse_std=eta*sigma0; %noise perterbation used in flipping function E

        n = sigma0*randn(1,N)+sigma0*randn(1,N)*1i; %complex noise
        y = c + n; % codeword+noise, channel information

        for j=1:N
            distance=abs(qam8-y(j));
            [~,min_idx]= min(distance);
            hard_d_cmplx(j) = qam8(min_idx); % hard decision in complex
            hard_d_gf8(j) = gf8(min_idx);
        end
        
        errors = hard_d_cmplx ~= c; 
        n_errors_hard =sum(errors); %total number of symbol errors after hard decision made
              
        errors_uncoded_bit = zeros(1, K);

        %bit error calculation in hard decision
        for e = 1 : K          
                if hard_d_gf8(e)~=code_seq(e) 
                    s1 = dec2bin(code_seq(e),p);
                    s2 = dec2bin(hard_d_gf8(e),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_uncoded_bit(e) = errors_uncoded_bit(e) + num_diff_bit;                
                end
        end
        
        un_bit_error = sum(errors_uncoded_bit); % no of bit error in each frame
        BE_unCoded(i)=BE_unCoded(i)+un_bit_error; % total no. of bit error in total frame generated in current Eb/No
        
        %calling the decoing function    
        [seqgf,failed,l]= decodeMultivote(code_seq,hard_d_cmplx, hard_d_gf8,...
        qam8, gf8,y, h, N, M, T, w,add_mat,mul_mat,div_mat,...
        CN_lst, nse_std,qam_binary_map,flip_num,No);  

               
        iters_cnr(i)=iters_cnr(i)+l;

        %bit error calculation for decoded sequence 
    
        errors_coded_bit = zeros(1, K);
            for g = 1 : K            
                if seqgf(g)~=code_seq(g)
                    s1 = dec2bin(code_seq(g),p);
                    s2 = dec2bin(seqgf(g),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_coded_bit(g) = errors_coded_bit(g) + num_diff_bit;                
                end     
            end

        bit_error = sum(errors_coded_bit); % no. of bit error in each frame
        BE_Coded(i)=BE_Coded(i)+bit_error; % total no. of bit error in total frame generated in the current Eb/No

        %iters_cnr(i) = iters_cnr(i)+l;

        frame_error = any(seqgf ~= code_seq);

        if frame_error
        
            FE(i) = FE(i) + 1;
            fail_count = fail_count + 1;
        
            % Error-symbol locations after decoding
            err_pos = find(seqgf ~= code_seq);
            a = length(err_pos);   % number of erroneous variable nodes
        
            % Final syndrome
            final_syn = decod_prod(seqgf, h, CN_lst, mul_mat, add_mat);
            unsat_cn = find(final_syn ~= 0);
            b = length(unsat_cn);  % number of unsatisfied check nodes


            % True and decoded GF(8) symbols at erroneous VNs
            true_err_symbols = code_seq(err_pos);
            decoded_err_symbols = seqgf(err_pos);
            
            % Error symbols: e_v = decoded_v - transmitted_v
            % In GF(2^p), subtraction is the same as addition
            error_symbols = zeros(1, length(err_pos));

            for ee = 1:length(err_pos)
            
                true_val = true_err_symbols(ee);
                dec_val  = decoded_err_symbols(ee);
            
                error_symbols(ee) = add_mat(dec_val + 1, true_val + 1);
            
            end

            % Check-node syndrome contribution details
            cn_details = struct([]);
            
            % CNs connected to at least one erroneous VN
            affected_cn = find(any(h(:, err_pos) ~= 0, 2));

            for cc = 1:length(affected_cn)
            
                cn = affected_cn(cc);
            
                % Erroneous VNs connected to this CN
                connected_error_vns = intersect(err_pos, CN_lst{cn});
            
                manual_syndrome = 0;
                contribution_values = zeros(1, length(connected_error_vns));
                edge_coefficients = zeros(1, length(connected_error_vns));
                connected_error_symbols = zeros(1, length(connected_error_vns));
            
                for vv = 1:length(connected_error_vns)
            
                    vn = connected_error_vns(vv);
            
                    % Position of this VN in err_pos
                    err_idx = find(err_pos == vn, 1);
            
                    % GF(8) error value
                    e_vn = error_symbols(err_idx);
            
                    % Edge coefficient h_{cn,vn}
                    h_cv = h(cn, vn);
            
                    % Contribution h_{cn,vn} * e_vn
                    contribution = mul_mat(h_cv + 1, e_vn + 1);
            
                    % Add contribution to syndrome in GF(8)
                    manual_syndrome = ...
                        add_mat(manual_syndrome + 1, contribution + 1);
            
                    edge_coefficients(vv) = h_cv;
                    connected_error_symbols(vv) = e_vn;
                    contribution_values(vv) = contribution;
            
                end
            
                cn_details(cc).cn = cn;
                cn_details(cc).error_vns = connected_error_vns;
                cn_details(cc).edge_coefficients = edge_coefficients;
                cn_details(cc).error_symbols = connected_error_symbols;
                cn_details(cc).contributions = contribution_values;
                cn_details(cc).manual_syndrome = manual_syndrome;
                cn_details(cc).decoder_syndrome = final_syn(cn);
            
            end

            

            % Store failed-frame information
            fail_log(fail_count).EbNo = Eb_No_db(i);
            fail_log(fail_count).frame = genFrame(i);
            fail_log(fail_count).iteration = l;
            fail_log(fail_count).initial_err_pos = find(hard_d_gf8 ~= code_seq);
            fail_log(fail_count).n_initial_errors = length(find(hard_d_gf8 ~= code_seq));
            fail_log(fail_count).err_pos = err_pos;
            fail_log(fail_count).unsat_cn = unsat_cn;
            fail_log(fail_count).a = a;
            fail_log(fail_count).b = b;
            fail_log(fail_count).seqgf = seqgf;
            fail_log(fail_count).code_seq = code_seq;
 
            % Additional symbol-level information
            fail_log(fail_count).true_err_symbols = true_err_symbols;
            fail_log(fail_count).decoded_err_symbols = decoded_err_symbols;
            fail_log(fail_count).error_symbols = error_symbols;
            fail_log(fail_count).final_syn = final_syn;

            fail_log(fail_count).affected_cn = affected_cn;
            fail_log(fail_count).cn_details = cn_details;

            
        end      
     end
     fprintf('Eb/No = %.1f dB: BER (coded) = %.6e, BER (uncoded) = %.6e\n', ...
        Eb_No_db(i), ...
        BE_Coded(i) / (genFrame(i) * K * p), ...
        BE_unCoded(i) / (genFrame(i) * K * p));
end

fprintf('\nDominant failure / trapping-set candidates:\n');

ab_pairs = zeros(length(fail_log),2);

for k = 1:length(fail_log)
    ab_pairs(k,:) = [fail_log(k).a, fail_log(k).b];
end

unique_pairs = unique(ab_pairs, 'rows');

for u = 1:size(unique_pairs,1)
    a = unique_pairs(u,1);
    b = unique_pairs(u,2);

    count = sum(ab_pairs(:,1)==a & ab_pairs(:,2)==b);

    fprintf('(%d,%d) occurred %d times\n', a, b, count);
end


pair_list = [];

for k = 1:length(fail_log)
    if fail_log(k).a == 2 && fail_log(k).b == 1
        pair_list = [pair_list; sort(fail_log(k).err_pos)];
    end
end

[unique_pairs, ~, idx] = unique(pair_list, 'rows');
counts = accumarray(idx, 1);

[counts_sorted, order] = sort(counts, 'descend');
unique_pairs_sorted = unique_pairs(order,:);

fprintf('\nMost frequent VN pairs for (2,1):\n');
for r = 1:min(10,length(counts_sorted))
    fprintf('VN pair (%d,%d) occurred %d times\n', ...
        unique_pairs_sorted(r,1), unique_pairs_sorted(r,2), counts_sorted(r));
end

vn1 = 26; %22;
vn2 = 31; %37;

cn_vn1 = find(h(:,vn1) ~= 0);
cn_vn2 = find(h(:,vn2) ~= 0);

fprintf('\nVN %d connected CNs: ', vn1);
fprintf('%d ', cn_vn1);

fprintf('\nVN %d connected CNs: ', vn2);
fprintf('%d ', cn_vn2);

shared_cn = intersect(cn_vn1, cn_vn2);

fprintf('\nShared CNs: ');
fprintf('%d ', shared_cn);
fprintf('\n');

%finding odd degree CN
T = [26 31];   % dominant VN pair

% Submatrix induced by these VNs
H_T = h(:, T);

% Induced degree of each CN within this VN subset
cn_deg = sum(H_T ~= 0, 2);

% CNs participating in the induced subgraph
induced_cn = find(cn_deg > 0);

% Odd-degree and even-degree CNs in the induced subgraph
odd_cn  = find(mod(cn_deg, 2) == 1);
even_cn = find(cn_deg > 0 & mod(cn_deg, 2) == 0);

fprintf('\nInduced CNs for VN pair (%d,%d): ', T(1), T(2));
fprintf('%d ', induced_cn);

fprintf('\nOdd-degree CNs in induced subgraph: ');
fprintf('%d ', odd_cn);

fprintf('\nEven-degree CNs in induced subgraph: ');
fprintf('%d ', even_cn);

fprintf('\nNumber of odd-degree CNs = %d\n', length(odd_cn));
fprintf('Number of even-degree CNs = %d\n', length(even_cn));

BERunCoded= BE_unCoded ./(genFrame * K *p); %bit error rate for uncoded

BERCoded = BE_Coded ./ (genFrame * K *p); %bit error rate for coded


 figure;

 semilogy(Eb_No_db, BERCoded, 'gx-', 'LineWidth', 1.2);
 grid on;
 hold on;

 semilogy(Eb_No_db, BERunCoded, 'rx-','LineWidth', 1.2);

 
 hold off;

 ylim([10e-8 10e-1]);
 xlim([2 20]);
 xlabel('E_b/N_0 (dB)', 'FontSize', 14);
 ylabel('BER', 'FontSize', 14);
 title('BER curve LDPC code (45,2,3) over GF(8)');
 hold off;
 legend( 'NG-CNV-SF', 'Uncoded','Location', 'northeast');







