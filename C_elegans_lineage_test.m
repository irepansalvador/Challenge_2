% simulate the deletions using parameters from the experimental data
% the data is based on the "fast" target in the embryo
% for the explanation see the notebook "Fast_sequence_Embryo_Reader.ipynb"
%-------------------------------------------------------------------------
% set the random seed to 123 so it can be reproducible
rng(123);

% Define the MAIN parameters
%-------------------------------------------------------------------------
% Number of targets
targets = 200;
% Mutation rate per minute
Lambda = 0.0002;
Norm_time = 1.88; % this is a normalising factor to get the real mins
%NUmber of sims
Nsims = 1;
% the number of states in the analysis. PAUP alows max 64 states.
states = 30;
% max_time
max_time = 3065;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for r = 1:Nsims                      % TO HAVE MULTIPLE SAMPLES
    %r = 1;
    disp("Rep=" + r);

    %%%%%%% CREATE A TABLE WITH THE PROBABILITY OF MUTATION IN EACH TARGET,
    %%%%%%% WHICH WILL FOLLOW A RANDOM GAMMA DISTRIBUTION
    rvs = gamrnd(0.1,2,[1 states]); rvs = rvs/sum(rvs);
    Nmer_prob = sort(rvs);
    Nmer_prob = cumsum(Nmer_prob);  
    % ---------------------------------------------------------------------
    %the array is a 2D array with the targets (within the cell) as columns and
    % each row is a daughter cell (matrix changes in each iteration)
    % ----------Initialize some variables ------------------------------------
    sequences = struct('Header',{},'Sequence',{});
    A1 = zeros(1,targets); % for the current cell division
    A2 = [];               % to store the cells of the next gen
    A3 = [];               % to store the terminal cells

    C_eleg = jsondecode(fileread('json-celllineage_MOD.js'));
    names = 'did';

    % this is the initial node that is goind to be analysed
    node = "C_eleg";

    found = 1;
    Ndiv = 1;
    daughters = node;
    N = 0;
    % make a recursive function that traverses the tree cell division by
    % cell division and only stops until there are no more cell daughters
    while found > 0
      f = 0;  % to store the total number of children in a givel cell div
      disp("    Division=" + Ndiv);
      Gdaughters = [];
      % -------- FOR EACH OF THE DAUGHTERS
      for d = 1:length(daughters)
        % print the node analysed
        %disp(daughters(d));
        % retrieves the data asociated with the node
        field = fieldnames(eval(daughters(d)))';
        % check if the node has the "children" field
        x = length(find(strcmp(field,'children')));
        % get the recorder for the cell
        Atemp = A1(d,:);
        % if the node has children  --------------------------------
        if x == 1
          f = f+1;
          % get number of children
          n_k = length(eval(daughters(d) + ".children"));
%          disp(daughters(d) + " has " + n_k +  " children" );
          for k = 1:n_k
            %N = N+1;
            % push every child to a new matrix
            if isa(eval(strcat(daughters(d),".children")),'cell')
              Gdx = strcat(daughters(d),".children{"+k+"}");
            else
              Gdx = strcat(daughters(d),".children("+k+")");
            end
            id = eval(strcat(Gdx, ".did") );
            t1  = round(eval(strcat(Gdx, ".levelDistance"))/1.88);
            t0  = round(eval(strcat(Gdx, ".totalDistance"))/1.88) - t1;
            t1  = t0 + t1;
            %disp(id + " born in " + t0 + "-" + t1 + " mins" );

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ----- simulate accumulation of muts (internal cells)
            events = zeros(targets,1);
            for i = 1:size(events)
              a = poisson_fixed_events(Lambda,1);
              events(i,1) = a(2);
            end
            Nevents = size(events(events > t0 & events <= t1),1);
            % for each copied target, mutate randomly with a fixed probability
            % the site can mutate only once, into one of n possible states
            Aind = find(Atemp == 0);
            AR = 0;
            if size(Aind,2) > 0 %%%% only if there are unmutated targets!
              if Nevents >= size(Aind,2)
                Eind = Aind;
              end

              if Nevents < size(Aind,2)
                Eind = randsample(Aind,Nevents);
              end

              AR = rand(size(Eind));
              AR = arrayfun(@(z)sum(z <= Nmer_prob), AR);

              Atemp(Eind) = AR; %%% MODIFY HEREEE!! Aind should be the ref
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ----------------------------------------
            A2 = vertcat(A2,Atemp);

            Gdaughters = horzcat(Gdaughters,Gdx);
          end
          % if the node has NO children -----------------------------
        else
          N = N+1;
          term_id = eval(strcat(daughters(d), ".did"));
          term_t1=round(eval(strcat(daughters(d),".levelDistance"))/1.88);
          term_t0=round(eval(strcat(daughters(d),".totalDistance"))/1.88);
          term_t0=term_t0 - term_t1;
          %term_t1 = term_t0 + term_t1;
          term_t1 = max_time ;
          disp(term_id+"born in "+ term_t0+". Muts accum until-"+term_t1);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % ----- simulate accumulation of muts (terminal cells)-
          events = zeros(targets,1);
          for i = 1:size(events)
            a = poisson_fixed_events(Lambda,1);
            events(i,1) = a(2);
          end
          Nevents = size(events(events > term_t0 & events <= term_t1),1);
          % for each copied target, mutate randomly with a fixed probability
          % the site can mutate only once, into one of n possible states
          Aind = find(Atemp == 0);
          AR = 0;
          if size(Aind,2) > 0 %%%% only if there are unmutated targets!
            if Nevents >= size(Aind,2)
              Eind = Aind;
            end

            if Nevents < size(Aind,2)
              Eind = randsample(Aind,Nevents);
            end

            AR = rand(size(Eind));
            AR = arrayfun(@(z)sum(z <= Nmer_prob), AR);

            Atemp(Eind) = AR;     %%% MODIFY HEREEE!! Aind should be the ref
          end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          A3 = vertcat(A3,Atemp);

          sequences(N).Header = term_id;
        % sequences(N).Sequence = Atemp;
        end
      end
      daughters = Gdaughters;
      A1 = A2;
      A2 = [];
      Ndiv = Ndiv + 1;
      if f == 0
        found = 0;
      end
    end
    disp("Total N cells = " + N);

    % Generate output
    % ------------- IO --------------------------------------------
    % name of the simulation
    sim = [sprintf('%02d',targets),'_targets_'...
        sprintf('%02d',states),'_states_'...
        sprintf('%.04f',Lambda),'_Lambda_'...
        ];
    A1 = char(zeros(size(A3)));
    %define the output file name for the alignment (fasta)
    out = ['C_elegans_',sim, sprintf('%04d', r),'rep','.fas'];
    freq_file = [sim, sprintf('%04d', r),'rep_C_elegans_freqs','.txt'];
    %---------- produce the output for each division!!
    %%%% Converts to ASCII characters, including 0-9, A-Z, a-z
    %%%% Avoids troublesome characters as ?, \, etc..
    % 26 most common (a-z)
    A1(A3< 27 & A3 >0 )= char(A3(A3< 27 & A3 >0) + 64);
    % Rare (A-Z)
    A1(A3>=27 & A3 <=30) = char(A3(A3>=27 & A3 <=30) + 70);
    %unmutated 
    A1(A3==0)= char(A3(A3==0) + 48);

    %%%%%% EXPORT THE SIMULATED FREQS TO A FILE %%%%%%%%%%%%
    unqA = unique(A1(1:numel(A1)));
    Nmers = cellstr(unqA');
    countNmers=histc(A1(1:numel(A1)),unqA); 
    relFreq=countNmers/numel(A1);
    %transp(relFreq);
    T = table(Nmers, transpose(relFreq));

    % create the file
    writetable(T,['./simulations/',freq_file], ...
    'WriteVariableNames',false);            

    myfasta = fopen(['./simulations/',out],'w');
    fclose(myfasta);

    myfasta = fopen(['./simulations/',out],'a');
    for ii = 1:N
        sequences(ii).Sequence = A1(ii,:);
    end
    %%%% SORT THE CELLS IN THE OUTPUT SO IT IS EASY TO READ
    % [~,index] = sortrows({sequences.Header}.');
    % sequences = sequences(index); clear index

    for ii = 1:N
        % append to a file
        fprintf(myfasta,'%s\t',sequences(ii).Header);
        fprintf(myfasta,'%s\n',sequences(ii).Sequence);
    end
    fclose(myfasta);
end