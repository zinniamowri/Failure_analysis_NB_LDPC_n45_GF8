%% CREATE A SHORT FAILURE REPORT
clear;
clc;

matFile = 'failed_frame_log_15p0dB.mat';
reportFile = 'report_for_target_pair.txt';

if ~isfile(matFile)
    error('Cannot find MAT file: %s', matFile);
end

S = load(matFile);

fid = fopen(reportFile, 'w');

if fid == -1
    error('Could not create report file.');
end

cleanupObject = onCleanup(@() fclose(fid));

%% =========
% TITLE
% ==========

fprintf(fid, 'FAILED-FRAME ANALYSIS AT %.1f dB\n', S.Eb_No_db);
fprintf(fid, '%s\n\n', repmat('=', 1, 55));

%% ======================
% 1. SIMULATION SUMMARY
% =======================

fprintf(fid, '1. SIMULATION SUMMARY\n');
fprintf(fid, '%s\n', repmat('-', 1, 55));

observedFER = double(S.fail_count) / double(S.genFrame);
matchingCount = numel(S.matching_frames);
matchingPercentage = 100 * matchingCount / S.fail_count;

fprintf(fid, 'Generated frames          : %d\n', S.genFrame);
fprintf(fid, 'Failed frames             : %d\n', S.fail_count);
fprintf(fid, 'Observed FER              : %.3e\n', observedFER);
fprintf(fid, 'Target residual VN pair   : (%d,%d)\n', ...
    S.target_pair(1), S.target_pair(2));
fprintf(fid, 'Failures with target pair : %d of %d (%.2f%%)\n\n', ...
    matchingCount, S.fail_count, matchingPercentage);



%% ======================================
% 2. SYMBOL TRANSITIONS FOR TARGET PAIR
% =======================================

fprintf(fid, '3. SYMBOL TRANSITIONS FOR THE TARGET PAIR\n');
fprintf(fid, '%s\n', repmat('-', 1, 55));

transitions = double(S.unique_transitions);
counts = double(S.transition_counts(:));

[counts, sortIndex] = sort(counts, 'descend');
transitions = transitions(sortIndex,:);

fprintf(fid, ...
    '  VN %d transition   VN %d transition   Count   Percentage\n', ...
    S.target_pair(1), S.target_pair(2));

fprintf(fid, '  %s\n', repmat('-', 1, 49));

for k = 1:size(transitions,1)
    percentage = 100 * counts(k) / matchingCount;

    fprintf(fid, ...
        '       %d -> %-2d           %d -> %-2d        %3d      %6.2f%%\n', ...
        transitions(k,1), ...
        transitions(k,2), ...
        transitions(k,3), ...
        transitions(k,4), ...
        counts(k), ...
        percentage);
end



%% ===========================
% 4. CHECK-NODE CALCULATION
% ============================

fprintf(fid, '4. CHECK-NODE CALCULATION\n');
fprintf(fid, '%s\n', repmat('-', 1, 55));

% Find a target-pair frame having the dominant transition
representativeIndex = [];

for k = 1:numel(S.fail_log)
    entry = S.fail_log(k);

    samePair = isequal( ...
        sort(double(entry.err_pos(:).')), ...
        sort(double(S.target_pair(:).')));

    sameTransition = ...
        numel(entry.true_err_symbols) == 2 && ...
        numel(entry.decoded_err_symbols) == 2 && ...
        isequal(double(entry.true_err_symbols(:).'), ...
                transitions(1,[1 3])) && ...
        isequal(double(entry.decoded_err_symbols(:).'), ...
                transitions(1,[2 4]));

    if samePair && sameTransition
        representativeIndex = k;
        break;
    end
end

if isempty(representativeIndex)
    fprintf(fid, 'No representative frame was found.\n');
else
    entry = S.fail_log(representativeIndex);

    fprintf(fid, 'frame number: %d\n', entry.frame);
    fprintf(fid, 'Final erroneous VNs            : %s\n', ...
        vectorToString(entry.err_pos));
    fprintf(fid, 'Transmitted symbols            : %s\n', ...
        vectorToString(entry.true_err_symbols));
    fprintf(fid, 'Decoded symbols                : %s\n', ...
    vectorToString(entry.decoded_err_symbols));
    fprintf(fid, 'Unsatisfied CN                 : %s\n\n', ...
        vectorToString(entry.unsat_cn));

   

    for j = 1:numel(entry.cn_details)
        detail = entry.cn_details(j);

        fprintf(fid, 'CN %d:\n', detail.cn);
        fprintf(fid, '  Connected erroneous VNs : %s\n', ...
            vectorToString(detail.error_vns));
        fprintf(fid, '  Edge coefficients       : %s\n', ...
            vectorToString(detail.edge_coefficients));
        
        fprintf(fid, '  Syndrome                : %d\n', ...
            detail.manual_syndrome);

       
    end
end
%% ============================================================
% 5. FAILURE-PATTERN SUMMARY AND COMPLETE DECODED FRAMES
% =============================================================

fprintf(fid, '\n5. FAILURE-PATTERN SUMMARY\n');
fprintf(fid, '%s\n', repmat('-', 1, 75));

numberOfFailures = numel(S.fail_log);

% Define a failure pattern using:
% final erroneous VNs, true symbols, decoded symbols, and unsatisfied CNs
failureSignatures = strings(numberOfFailures,1);

for k = 1:numberOfFailures
    entry = S.fail_log(k);

    failureSignatures(k) = sprintf( ...
        'FinalVNs=%s|True=%s|Decoded=%s|UnsatCN=%s', ...
        vectorToString(entry.err_pos), ...
        vectorToString(entry.true_err_symbols), ...
        vectorToString(entry.decoded_err_symbols), ...
        vectorToString(entry.unsat_cn));
end

% Group identical failure patterns
[uniqueSignatures, ~, groupNumber] = unique( ...
    failureSignatures, 'stable');

patternCount = accumarray(groupNumber, 1);

fprintf(fid, ...
    '%7s %-12s %-12s %-12s %-10s\n', ...
    'Count', ...
    'Final VNs idx', ...
    'True sym', ...
    'Decoded sym', ...
    'Unsat CN');

fprintf(fid, '%s\n', repmat('-', 1, 75));

for g = 1:numel(uniqueSignatures)

    indices = find(groupNumber == g);
    representativeEntry = S.fail_log(indices(1));

    fprintf(fid, ...
        '%7d %-12s %-12s %-12s %-10s\n', ...
        patternCount(g), ...
        vectorToString(representativeEntry.err_pos), ...
        vectorToString(representativeEntry.true_err_symbols), ...
        vectorToString(representativeEntry.decoded_err_symbols), ...
        vectorToString(representativeEntry.unsat_cn));
end

%% ============================================================
% COMPLETE DECODED FRAME FOR EACH UNIQUE FAILURE PATTERN
% =============================================================

fprintf(fid, '\n');
fprintf(fid, '6. COMPLETE DECODED FRAME FOR EACH UNIQUE PATTERN\n');
fprintf(fid, '%s\n', repmat('-', 1, 75));

for g = 1:numel(uniqueSignatures)

    indices = find(groupNumber == g);
    representativeIndex = indices(1);
    entry = S.fail_log(representativeIndex);

    fprintf(fid, '\nPattern %d\n', g);
    fprintf(fid, '%s\n', repmat('-', 1, 40));

    fprintf(fid, 'Number of occurrences     : %d\n', patternCount(g));
    %fprintf(fid, 'Representative log index  : %d\n', representativeIndex);
    %fprintf(fid, 'Generated frame number    : %d\n', entry.frame);

    fprintf(fid, 'Final erroneous VNs       : %s\n', ...
        vectorToString(entry.err_pos));

    fprintf(fid, 'True symbols at errors    : %s\n', ...
        vectorToString(entry.true_err_symbols));

    fprintf(fid, 'Decoded symbols at errors : %s\n', ...
        vectorToString(entry.decoded_err_symbols));

    fprintf(fid, 'Unsatisfied CNs           : %s\n', ...
        vectorToString(entry.unsat_cn));

    % Print all generated frame numbers belonging to this pattern
    generatedFrames = zeros(1,numel(indices));

    for j = 1:numel(indices)
        generatedFrames(j) = S.fail_log(indices(j)).frame;
    end

    %fprintf(fid, 'Generated frame numbers   : %s\n', ...
      %vectorToString(generatedFrames));

    %% Print complete decoded frame
    fprintf(fid, '\nComplete decoded frame:\n');

    if isfield(entry, 'seqgf') && ~isempty(entry.seqgf)

        printLongVector(fid, entry.seqgf, 15);

    elseif isfield(entry, 'decoded_frame') && ...
            ~isempty(entry.decoded_frame)

        printLongVector(fid, entry.decoded_frame, 15);

    else
        fprintf(fid, ...
            ['Complete decoded frame was not stored in fail_log.\n', ...
             'Only decoded_err_symbols is available.\n']);
    end
end



%% ===============
% LOCAL FUNCTION
% ================

function output = vectorToString(value)

    if isempty(value)
        output = '[]';
        return;
    end

    value = double(value(:).');

    output = sprintf('%g,', value);
    output(end) = [];

    output = ['[', output, ']'];
end

function printLongVector(fid, values, valuesPerLine)

    values = double(values(:).');
    numberOfValues = numel(values);

    for k = 1:numberOfValues

        fprintf(fid, '%3d', values(k));

        if mod(k, valuesPerLine) == 0 || k == numberOfValues
            fprintf(fid, '\n');
        else
            fprintf(fid, ' ');
        end
    end
end
