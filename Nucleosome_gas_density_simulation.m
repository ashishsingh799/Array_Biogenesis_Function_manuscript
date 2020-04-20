function nucleosomegas_NatComm(fResidualHistones, sigma, numsimu, genelength_minmax, nDissociate, prctDissociate)
%
% Simulates statistical depositon of nucleosomes next to a well positioned +1 nucleosome.
% Nucleosomes are initially in a 'gas', ie unbound, and then settle one by
% one on DNA molecules that have a fixed +1 nucleosome. After an initial
% round, a fraction of nucleosomes is randomly selected and dissociated
% again, before thos enuucleosomes randomly settle on free stretches of DNA.
%
% Paramters:
% fResidualHistones: fraction of residual histones after histone depletion.
%                    0.6 refers to a reduction of histone content by 40% relative to wild-type
% sigma: refers to the sigma of a gaussian that simulates the well positioned +1 nucleosome
% numsimu: number of simulated DNA molecules
% genelength_minmax: is a vector [min max] defining the minimal and maximal lengths of the DNA moelcules in the simulation
% nDissociate : number of times, a certain percentage of nucleosomes is randomly dissociated and re-depositedprctDissociate)
% prctDissociate: percentage of nucleosomes that is randomly dissociated
%
% Example usage: nucleosomegas_NatComm(0.6, 25, 50000, [2000, 2500], 10, 20)
%
% (c) 2019 - 2020, Felix Mueller-Planitz

footprint = 146;    % nucleosome footprint in bp
wt_linker=19;       % linekr length in bp in WT yeast

stop_fractionbound = 0.03;  % when LESS than this fraction of nucs is in the gas ie not bound , quit simulation
stop_countUnsuccessful = 10000;   % after this many trials to place a nuc unsuccessfully, stop simu

% -- calculate linker length, avg occupancy
size_genome = 12000000;                             % doesnt matter what value. Cancels out. 12 Mbp are for S. cerevisiae
numnucs_wt = size_genome / (footprint+wt_linker);   % number of nucleosomes in WT. ~73000 for 12 Mbp genome

numnucs_hd = numnucs_wt * fResidualHistones;        % number of nucleosomes afte rhistone depletion
avg_occupancy = numnucs_hd * footprint / size_genome;

% -- make normal distribution to draw the position of the first nucleosome
pd = makedist('Normal', 'mu',2*sigma, 'sigma',sigma);
gaussian= floor(random(pd,numsimu,1));

% -- place the first nucleosome on each simulated DNA molecule
for ii = numsimu:-1:1
    dyad{ii}(1)= gaussian(ii);    % positzion the first nucleosome around a gaussian distribution
end

% -- make genes with predetermined sizes 'genelength_minmax' and calculate # of nucs in the gas
genelengths = genelength_minmax(1) + randi(diff(genelength_minmax), numsimu, 1);
max_nucs = ceil(avg_occupancy * sum(genelengths) / footprint);
nucs_in_gas = max_nucs - numsimu;   % number of not-yet-positioned nucs (the +1 is positioned already)


% --------------------------------------------------------------------
% simulate dyad positiones as a 1D gas for 'numsimu' number of genes
% --------------------------------------------------------------------
count = 0;
oldtime = cputime;
tic

for iDissociate = 1:nDissociate
    % ------ get rid of x% of all nucs_:
    if nucs_in_gas < max_nucs - numsimu     % this condition prevents it that during the first round nucs are dissociated already even before they are placed
        % but make sure that the +1 nucleosome is never dissociated
        nucs_placed = max_nucs - nucs_in_gas;
        for ii = 1:ceil(nucs_placed/100*prctDissociate)
            [nucs_in_gas, dyad] = eliminateNuc(dyad, nucs_in_gas, gaussian);
        end
    end
    
    % ------ settle nucleosomes in gas onto DNA:
    while nucs_in_gas > max_nucs * stop_fractionbound && count < stop_countUnsuccessful
        
        gene = randi(length(genelengths));  % pick random gene
        newdyad = randi(genelengths(gene)); % make random dyad position
        
        tmp_dyads = sort( [dyad{gene} newdyad] );
        NRL = diff( tmp_dyads  );
        if all( NRL >= footprint )    % make sure that all dyads are > 146 bp apart
            dyad{gene}= tmp_dyads;
            nucs_in_gas = nucs_in_gas - 1;
            count = 0;
        else
            % placemnet was unsuccessful
            count = count+1;
        end
    end
    
    if cputime - oldtime > 1    % display infos
        oldtime = cputime;
        disp(['Dissociation round #' num2str(iDissociate) ' completed. Elapsedtime= ' num2str(toc)]);
    end
end

% -- dyad composite plot
dyad_all = [dyad{:}];
dyad_binned=[];
dyad_binned(1,:) = 1:max(dyad_all); % make bin size = 1 bp
[dyad_binned(2,:)]= histc(dyad_all, dyad_binned(1,:));

plot(dyad_binned(1,:), smooth(  dyad_binned(2,:), 5  ), 'b-')

xlimits= [0, min(genelength_minmax)];
xlim(xlimits)



% -------------- subfunction -------------- 
% gets rid of a certain percentage of nucleosomes on random
% positions, but keep all +1 nucleosomes
function [nucs_in_gas, dyad] = eliminateNuc(dyad, nucs_in_gas, gaussian)
tmp_nnucs = cellfun('length', dyad);
unluckynuc = randi(sum(tmp_nnucs));   % pick random nuc to eliminate

cumsum_nucs = cumsum(tmp_nnucs);
idx =unluckynuc <= cumsum_nucs;
eliminate_gene = find(cumsum(idx)==1);
if eliminate_gene == 1
    eliminate_nuc = unluckynuc;
else
    eliminate_nuc = unluckynuc - cumsum_nucs(eliminate_gene-1);
end
if dyad{eliminate_gene}(eliminate_nuc) ~= gaussian(eliminate_gene)  % only release nucs if not at +1
    dyad{eliminate_gene}(eliminate_nuc)= [];
    nucs_in_gas = nucs_in_gas + 1;
end
