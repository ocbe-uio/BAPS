function baps
  version = [6 0 0 9001];
  versionStr = sprintf('%d.%d.%d.%d', version(1), version(2), version(3), version(4));
  disp(['Welcome to BAPS ' versionStr]);
  while true
    prompt = ['Please select the function you want to run ' ...
      '(1: greedyMix, 2: greedyPopMix): '];
    choice = input(prompt, 's');
    if isempty(choice)
      disp('Exiting BAPS');
      return;
    end
    switch choice
      case '1'
        disp('Clustering of individuals');
        greedyMix(-1);
        break;
      case '2'
        disp('Clustering of groups of individuals');
        greedyPopMix;
        break;
      otherwise
        disp('Invalid choice. Please try again.');
    end
  end
end
