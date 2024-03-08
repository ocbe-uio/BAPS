function baps(file, file_type, analysis)
  % Welcome message
  ver = [6 0 0 9001];
  versionStr = sprintf('%d.%d.%d.%d', ver(1), ver(2), ver(3), ver(4));
  disp(['Welcome to BAPS ' versionStr]);

  % Requesting analysis if not provided
  if nargin < 2 || isempty(analysis)
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

  % Dispatching analysis
  switch analysis
    case 'greedyMix'
      disp('Clustering of individuals');
      greedyMix(file, file_type);
    case 'greedyPopMix'
      disp('Clustering of groups of individuals');
      greedyPopMix;
    otherwise
      disp('Invalid choice. Please try again.');
  end
end
