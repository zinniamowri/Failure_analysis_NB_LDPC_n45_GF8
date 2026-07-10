%% SIMPLE ERROR-SYMBOL AND VN-PAIR ANALYSIS
clear;
clc;

matFile = 'failed_frame_log_16p0dB.mat';

if ~isfile(matFile)
    error('Cannot find MAT file: %s', matFile);
end

S = load(matFile);

if ~isfield(S, 'fail_log') || isempty(S.fail_log)
    error('The MAT file does not contain a valid fail_log.');
end

if ~isfield(S, 'Eb_No_db')
    error('The MAT file does not contain Eb_No_db.');
end

if ~isfield(S, 'genFrame')
    error('The MAT file does not contain genFrame.');
end

% Create a filename such as:
% error_symbol_pair_report_16p0dB.txt
snrText = strrep(sprintf('%.1f', S.Eb_No_db), '.', 'p');

reportFile = sprintf( ...
    'error_symbol_pair_report_%sdB.txt', snrText);

fid = fopen(reportFile, 'w');

if fid == -1
    error('Could not create report file: %s', reportFile);
end

cleanupObject = onCleanup(@() fclose(fid));

numberOfFailures = numel(S.fail_log);
observedFER = double(numberOfFailures) / double(S.genFrame);

%% =========================================
% REPORT TITLE AND SIMULATION INFORMATION
% ==========================================

fprintf(fid, 'ERROR-SYMBOL AND VN-PAIR ANALYSIS AT %.1f dB\n', ...
    S.Eb_No_db);

fprintf(fid, '%s\n\n', repmat('=', 1, 75));

fprintf(fid, 'SIMULATION INFORMATION\n');
fprintf(fid, '%s\n', repmat('-', 1, 75));

fprintf(fid, 'Eb/N0                     : %.1f dB\n', ...
    S.Eb_No_db);

fprintf(fid, 'Total generated frames    : %d\n', ...
    S.genFrame);

fprintf(fid, 'Total failed frames       : %d\n', ...
    numberOfFailures);

fprintf(fid, 'Observed frame error rate : %.6e\n\n', ...
    observedFER);

%% ==============================================
% 1. PRINT ERROR SYMBOLS FOR EVERY FAILED FRAME
% ===============================================

fprintf(fid, '1. ERROR SYMBOLS IN EVERY FAILED FRAME\n');
fprintf(fid, '%s\n', repmat('-', 1, 75));

fprintf(fid, ...
    '%10s %10s %15s %15s\n', ...
    'Frame', 'VN index', 'True symbol', 'Decoded symbol');

fprintf(fid, '%s\n', repmat('-', 1, 75));

for k = 1:numberOfFailures

    entry = S.fail_log(k);

    errorVNs = double(entry.err_pos(:).');
    trueSymbols = double(entry.true_err_symbols(:).');
    decodedSymbols = double(entry.decoded_err_symbols(:).');

    numberOfErrors = numel(errorVNs);

    if numel(trueSymbols) ~= numberOfErrors || ...
            numel(decodedSymbols) ~= numberOfErrors

        warning(['Frame %d has inconsistent numbers of VNs, true ', ...
                 'symbols, and decoded symbols.'], entry.frame);
        continue;
    end

    for j = 1:numberOfErrors

        fprintf(fid, ...
            '%10d %10d %15d %15d\n', ...
            entry.frame, ...
            errorVNs(j), ...
            trueSymbols(j), ...
            decodedSymbols(j));
    end
end

%% ============================
% 2. COLLECT ALL VN PAIRS
% =============================

% Each row will contain:
%
% column 1: VN 1
% column 2: VN 2
% column 3: true symbol at VN 1
% column 4: decoded symbol at VN 1
% column 5: true symbol at VN 2
% column 6: decoded symbol at VN 2
% column 7: generated frame number
% column 8: fail-log index

pairRecords = [];

for k = 1:numberOfFailures

    entry = S.fail_log(k);

    errorVNs = double(entry.err_pos(:).');
    trueSymbols = double(entry.true_err_symbols(:).');
    decodedSymbols = double(entry.decoded_err_symbols(:).');

    numberOfErrors = numel(errorVNs);

    if numberOfErrors < 2
        continue;
    end

    if numel(trueSymbols) ~= numberOfErrors || ...
            numel(decodedSymbols) ~= numberOfErrors
        continue;
    end

    % Sort the VNs so that, for example, (26,31) and (31,26)
    % are treated as the same pair.
    [errorVNs, sortIndex] = sort(errorVNs);

    trueSymbols = trueSymbols(sortIndex);
    decodedSymbols = decodedSymbols(sortIndex);

    % Generate all two-VN combinations in this failed frame.
    pairCombinations = nchoosek(1:numberOfErrors, 2);

    for j = 1:size(pairCombinations,1)

        firstIndex = pairCombinations(j,1);
        secondIndex = pairCombinations(j,2);

        newRecord = [ ...
            errorVNs(firstIndex), ...
            errorVNs(secondIndex), ...
            trueSymbols(firstIndex), ...
            decodedSymbols(firstIndex), ...
            trueSymbols(secondIndex), ...
            decodedSymbols(secondIndex), ...
            double(entry.frame), ...
            k];

        pairRecords = [pairRecords; newRecord]; 
    end
end

%% ========================================
% 3. COUNT HOW OFTEN EACH VN PAIR OCCURS
% =========================================

fprintf(fid, '\n\n2. VN-PAIR OCCURRENCE SUMMARY\n');
fprintf(fid, '%s\n', repmat('-', 1, 75));

if isempty(pairRecords)

    fprintf(fid, 'No failed frame contained two or more erroneous VNs.\n');

else

    allPairs = pairRecords(:,1:2);

    [uniquePairs, ~, pairGroup] = unique(allPairs, 'rows');

    pairCounts = accumarray(pairGroup, 1);

    % Number of distinct failed frames containing each pair
    pairFrameCounts = zeros(size(pairCounts));

    for g = 1:size(uniquePairs,1)

        groupRows = pairGroup == g;
        pairFrameCounts(g) = numel(unique(pairRecords(groupRows,7)));
    end

    % Sort pairs from most frequent to least frequent
    [pairCounts, sortIndex] = sort(pairCounts, 'descend');

    uniquePairs = uniquePairs(sortIndex,:);
    pairFrameCounts = pairFrameCounts(sortIndex);

    fprintf(fid, ...
        '%8s %10s %10s %15s %15s\n', ...
        'Rank', 'VN 1', 'VN 2', ...
        'Pair records', 'Failed frames');

    fprintf(fid, '%s\n', repmat('-', 1, 75));

    for g = 1:size(uniquePairs,1)

        fprintf(fid, ...
            '%8d %10d %10d %15d %15d\n', ...
            g, ...
            uniquePairs(g,1), ...
            uniquePairs(g,2), ...
            pairCounts(g), ...
            pairFrameCounts(g));
    end
end

%% =============================================
% 4. TRUE AND DECODED SYMBOLS FOR EACH VN PAIR
% ==============================================

fprintf(fid, '\n\n3. SYMBOL TRANSITIONS FOR EACH VN PAIR\n');
fprintf(fid, '%s\n', repmat('-', 1, 75));

if ~isempty(pairRecords)

    for g = 1:size(uniquePairs,1)

        vn1 = uniquePairs(g,1);
        vn2 = uniquePairs(g,2);

        pairRows = ...
            pairRecords(:,1) == vn1 & ...
            pairRecords(:,2) == vn2;

        currentRecords = pairRecords(pairRows,:);

        % Symbol pattern:
        %
        % true VN1, decoded VN1, true VN2, decoded VN2
        symbolPatterns = currentRecords(:,3:6);

        [uniqueSymbolPatterns, ~, symbolGroup] = ...
            unique(symbolPatterns, 'rows');

        symbolCounts = accumarray(symbolGroup, 1);

        [symbolCounts, symbolSortIndex] = ...
            sort(symbolCounts, 'descend');

        uniqueSymbolPatterns = ...
            uniqueSymbolPatterns(symbolSortIndex,:);

        fprintf(fid, '\nVN pair (%d,%d)\n', vn1, vn2);
        fprintf(fid, 'Pair occurrence count : %d\n', pairCounts(g));
        fprintf(fid, ...
            'Failed-frame count    : %d\n\n', pairFrameCounts(g));

        fprintf(fid, ...
            '  VN %d         VN %d         Count   Percentage\n', ...
            vn1, vn2);

        fprintf(fid, ...
            '  True->Decoded True->Decoded\n');

        fprintf(fid, '  %s\n', repmat('-', 1, 55));

        totalForPair = sum(symbolCounts);

        for j = 1:size(uniqueSymbolPatterns,1)

            percentage = ...
                100 * symbolCounts(j) / totalForPair;

            fprintf(fid, ...
                '      %d -> %-2d       %d -> %-2d      %5d      %6.2f%%\n', ...
                uniqueSymbolPatterns(j,1), ...
                uniqueSymbolPatterns(j,2), ...
                uniqueSymbolPatterns(j,3), ...
                uniqueSymbolPatterns(j,4), ...
                symbolCounts(j), ...
                percentage);
        end
    end
end

%% =============================
% 5. MOST CONTRIBUTING PAIR
% ==============================

fprintf(fid, '\n\n4. DOMINANT VN PAIR\n');
fprintf(fid, '%s\n', repmat('-', 1, 75));

if ~isempty(pairRecords)

    dominantPair = uniquePairs(1,:);

    fprintf(fid, 'Most frequently observed VN pair : (%d,%d)\n', ...
        dominantPair(1), dominantPair(2));

    fprintf(fid, 'Number of pair occurrences       : %d\n', ...
        pairCounts(1));

    fprintf(fid, 'Number of failed frames          : %d of %d\n', ...
        pairFrameCounts(1), numberOfFailures);

    fprintf(fid, 'Fraction of failed frames         : %.2f%%\n', ...
        100 * pairFrameCounts(1) / numberOfFailures);
end

fprintf('\nAnalysis completed.\n');
fprintf('Report written to: %s\n', reportFile);