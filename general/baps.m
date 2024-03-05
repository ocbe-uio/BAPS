function baps
  disp('Welcome to BAPS');
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
        greedyMix(-1);
        break;
      case '2'
        greedyPopMix;
        break;
      otherwise
        disp('Invalid choice. Please try again.');
    end
  end
end
