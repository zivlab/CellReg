function display_progress_bar(percent,terminate_bar)

% This function creates a text progress bar. It should be called with a
% STRING argument to initialize and terminate. Otherwise the number correspoding
% to progress in % should be supplied.
% INPUTS:   percent   Either: Text string to initialize or terminate
%                       Percentage number to show progress
% terminate_bar: true if a previous bars should be terminated first
% and false if not

% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m

% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version

% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/

%% Initialization
persistent strCR;           %   Carriage return pesistent variable

% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar

%% Main

if terminate_bar
    strCR = [];    
else
    if isempty(strCR) && ~ischar(percent),
        % Progress bar must be initialized with a string
        error('The text progress must be initialized with a string');
    elseif isempty(strCR) && ischar(percent),
        % Progress bar - initialization
        fprintf('%s',percent);
        strCR = -1;
    elseif ~isempty(strCR) && ischar(percent),
        % Progress bar  - termination
        strCR = [];
        fprintf([percent '\n']);
    elseif isnumeric(percent)
        % Progress bar - normal progress
        c = floor(percent);
        percentageOut = [num2str(percent) '%%'];
        percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
        nDots = floor(c/100*strDotsMaximum);
        dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
        strOut = [percentageOut dotOut];
        
        % Print it on the screen
        if strCR == -1,
            % Don't do carriage return during first run
            fprintf(strOut);
        else
            % Do it during all the other runs
            fprintf([strCR strOut]);
        end
        
        % Update carriage return
        strCR = repmat('\b',1,length(strOut)-1);
        
    else
        % Any other unexpected input
        error('Unsupported argument type');
    end
end

end