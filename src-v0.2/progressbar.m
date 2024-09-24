function progressbar(c)
    % textprogressbar() displays a text-based progress bar with an estimated time.
    % 
    % Use:
    % textprogressbar('Start Message'); % Initialize the progress bar
    % textprogressbar(P); % Update the progress bar at P (0 to 100)
    % textprogressbar(100); % Finalize the progress bar
    % 
    % Author: Selim Romero
    % 
    
    persistent startTime prevC avgElapsedTimes blockCount;
    barLength = 25; % Length of the progress bar

    if ischar(c)
        % Initialization with custom message
        fprintf('%s\n', c);
        startTime = tic; % Start timing
        prevC = 0; % Reset previous progress
        avgElapsedTimes = []; % Initialize array for block times
        blockCount = 0; % Reset block count
        return;
    elseif isnumeric(c)
        if isempty(startTime)
            error('Progress bar not initialized. Call with a string first.');
        end

        elapsedTime = toc(startTime); % Time since start
        progress = max(0, min(100, c)); % Ensure within bounds [0, 100]

        if progress == 100
            fprintf('100%% [%-*s] Done!\n', barLength, repmat('.', 1, barLength));
            prevC = 0;
            return;
        end
        
        % Only update if progress has changed significantly
        if abs(progress - prevC) < 1
            return;
        end

        % Record time taken for this block
        blockCount = blockCount + 1;
        avgElapsedTimes(end + 1) = elapsedTime; % Append current elapsed time

        % Calculate weighted average elapsed time
        if blockCount > 1
            weights = linspace(1, 0.5, length(avgElapsedTimes)); % Decreasing weights
            weightedAvgTime = sum(avgElapsedTimes .* weights) / sum(weights);
        else
            weightedAvgTime = elapsedTime; % Initial case
        end

        remainingTime = weightedAvgTime * (100 - progress) / progress; % Smoothed ETA
        timeStr = sec2timestr(remainingTime); % Format remaining time
        
        % Clear the line and update progress bar
        fprintf('\r%.0f%% [%-*s] ETA: %s', progress, barLength, repmat('.', 1, round(barLength * progress / 100)), timeStr);
        
        prevC = progress; % Update previous progress
    else
        error('Unsupported input type.');
    end
end

function timeStr = sec2timestr(sec)
    % Converts seconds into a readable format (e.g., hr, min, sec)
    h = floor(sec / 3600);
    m = floor(mod(sec, 3600) / 60);
    s = mod(sec, 60);

    if h > 0
        timeStr = sprintf('%dh %dm', h, m);
    elseif m > 0
        timeStr = sprintf('%dm %ds', m, floor(s));
    else
        timeStr = sprintf('%ds', floor(s));
    end
end

