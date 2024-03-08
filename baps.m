function baps(file, file_type, analysis, partitionCompare)
  % Adding functions from the current directory and its subdirectories
  addpath(genpath(cd));

  % Welcome message
  ver = [6 0 0 9001];
  versionStr = sprintf('%d.%d.%d.%d', ver(1), ver(2), ver(3), ver(4));
  disp(['Welcome to BAPS ' versionStr]);

  % Requesting analysis if not provided
  if nargin < 3 || isempty(analysis)
    prompt = ['Please select the function you want to run ' ...
      '(1: greedyMix, 2: greedyPopMix): '];
    analysis = input(prompt, 's');
    switch analysis
      case '1'
        analysis = 'greedyMix';
      case '2'
        analysis = 'greedyPopMix';
      otherwise
        disp('Invalid choice. Please try again.');
    end
    if isempty(analysis)
      disp('Exiting BAPS');
      return;
    end
  end

  % Processing boolean input
  partitionCompare = process_boolean_input(4);

  % Dispatching analysis
  switch analysis
    case 'greedyMix'
      disp('Clustering of individuals');
      greedyMix(file, file_type, partitionCompare);
    case 'greedyPopMix'
      disp('Clustering of groups of individuals');
      greedyPopMix;
    otherwise
      disp('Invalid choice. Please try again.');
  end
end
