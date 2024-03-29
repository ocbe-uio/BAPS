function [data_clique, data_separator, noalle_clique, noalle_separator] = transform2(data, component_mat, linkage_model)
% Filename: transform.m
% [data_clique,data_clique,noalle_clique, noalle_separator] = tranform2(data, component_mat, linkage_model)
%
% Description:
% Function for tranforming the profile data into clique data and separator data.
% For the linear graphical models, in each component, ncliques + nseparators = 2*nnodes - 3
% holds for each block(component in the graph).

% Author: Jing Tang
% Modified date: 18/10/2005
% The function for treating missing data is added.

% Input:
% data: the allelic profile data. Last col is the index vector.
% component_mat: the component matrix.
% isdependent:
% islinear:

% Output:
% data_clique: [nindividuals x (ncliques+nsingletons)]
% data_separator:  [nindividual x nseparators]


if strcmp(linkage_model,'linear')
    nindividuals = size(data,1);
    nloci = size(data,2)-1;
    % if ~all(diag(adj_mat)), adj_mat = adj_mat + diag(ones(nloci,1));, end;
    % [p,p,r]=dmperm(adj_mat);
    
    % UNSOLVED --- Need to find a better algorithm to count the numbers of
    % separators ,cliques and singletons.
    ncomponents = size(component_mat,1);
    cardinality = sum(component_mat>0,2)';
    ncliques = sum(cardinality-1);
    nsingletons = sum(cardinality==1);
    nseparators = sum(cardinality(find(cardinality>1))-2);
    
    data_clique = zeros(nindividuals, ncliques+nsingletons);
    noalle_clique = zeros(ncliques+nsingletons,1);
    data_separator = zeros(nindividuals, nseparators);
    noalle_separator = zeros(nseparators,1);
    % data_separator = []; % used for vectorization
    k = 1;
    for i = 1:ncomponents
        if (cardinality(i)==1) % singleton
            data_clique(:,k) = data(:,component_mat(i,1));
            temp = unique(data(:,component_mat(i,1)));
            if temp(1) == 0 % missing data found
                noalle_clique(k) = size(temp,1) -1;
            else
                noalle_clique(k) = size(temp,1);
            end
            k = k+1;
        else
            for j = 1:cardinality(i)-1
                
                % dealing with missing data
                % nonmissing = find(all(data(:,[component_mat(i,j) component_mat(i,j+1)]),2));
                % [b,m,n] = unique(data(nonmissing,[component_mat(i,j) component_mat(i,j+1)]),'rows');
                % data_clique(nonmissing,k) = n;
                % noalle_clique(k) = size(unique(data(nonmissing,component_mat(i,j))),1)*...
                %    size(unique(data(nonmissing,component_mat(i,j+1))),1);
                
                [b,m,n] = unique(data(:,[component_mat(i,j) component_mat(i,j+1)]),'rows');
                data_clique(:,k) = n;
                
                noalle_clique(k) = size(unique(data(:,component_mat(i,j))),1)*...
                    size(unique(data(:,component_mat(i,j+1))),1);
                k = k+1;
            end
            
        end
    end
    
    k = 1;
    for i = 1:ncomponents
        if (cardinality(i)>2)
            for j = 2:cardinality(i)-1
                data_separator(:,k) = data(:,component_mat(i,j));
                k = k+1;
            end
        end
    end
    if (k-1)~=nseparators
        error('ERROR in transform.m');
        return
    end
    for i = 1:nseparators
        % nonmissing = find(all(data_separator(:,i),2));
        % noalle_separator(i) = length(unique(data_separator(nonmissing,i)));
        noalle_separator(i) = length(unique(data_separator(:,i)));
    end
    
elseif strcmp(linkage_model,'codon')
    nindividuals = size(data,1);
    nloci = size(data,2)-1;
    ncomponents = size(component_mat,1);
    cardinality = sum(component_mat>0,2)';
    
    ncliques = sum(cardinality(find(cardinality>2))-2);
    nsingletons = sum(cardinality==1 | cardinality ==2);
    nseparators = sum(cardinality(find(cardinality>2))-3);
    
    data_clique = zeros(nindividuals,ncliques+nsingletons);
    noalle_clique = zeros(ncliques+nsingletons,1);
    data_seperator = zeros(nindividuals,nseparators);
    noalle_separator = zeros(nseparators,1);
    k = 1;
    for i = 1:ncomponents
        if (cardinality(i)==1) % singleton
            data_clique(:,k) = data(:,component_mat(i,1));
            temp = unique(data(:,component_mat(i,1)));
            if temp(1) == 0 % missing data found
                noalle_clique(k) = size(temp,1) -1;
            else
                noalle_clique(k) = size(temp,1);
            end
            k = k+1;
        elseif (cardinality(i)==2) % transform to linear linkage
            % dealing with missing data
            % nonmissing = find(all(data(:,[component_mat(i,1) component_mat(i,2)]),2));
            % [b,m,n] = unique(data(nonmissing,[component_mat(i,1) component_mat(i,2)]),'rows');
            
            [b,m,n] = unique(data(:,[component_mat(i,1) component_mat(i,2)]),'rows');
            data_clique(:,k) = n;
            % data_clique(nonmissing,k) = n;
            
            % noalle_clique(k) = size(unique(data(nonmissing,component_mat(i,1))),1)*...
            % size(unique(data(nonmissing,component_mat(i,2))),1);
            noalle_clique(k) = size(unique(data(:,component_mat(i,1))),1)*...
                size(unique(data(:,component_mat(i,2))),1);
            k = k+1;
        else
            for j = 1:cardinality(i)-2
                % dealing with missing data
                % nonmissing = find(all(data(:,[component_mat(i,j) component_mat(i,j+1) component_mat(i,j+2)]),2));
                % [b,m,n] = unique(data(nonmissing,[component_mat(i,j) component_mat(i,j+1) component_mat(i,j+2)]),'rows');
                % data_clique(nonmissing,k) = n;
                
                [b,m,n] = unique(data(:,[component_mat(i,j) component_mat(i,j+1) component_mat(i,j+2)]),'rows');
                data_clique(:,k) = n;
                
                % noalle_clique(k) = size(unique(data(nonmissing,component_mat(i,j))),1)*...
                %    size(unique(data(nonmissing,component_mat(i,j+1))),1)*size(unique(data(nonmissing,component_mat(i,j+2))),1);
                noalle_clique(k) = size(unique(data(:,component_mat(i,j))),1)*...
                    size(unique(data(:,component_mat(i,j+1))),1)*size(unique(data(:,component_mat(i,j+2))),1);
                k = k+1;     
            end
        end
    end
    
    k = 1;
    for i = 1:ncomponents
        if (cardinality(i)>2)
            for j = 2:cardinality(i)-2
                % dealing with missing data
                % nonmissing = find(all(data(:,[component_mat(i,j) component_mat(i,j+1)]),2));
                % [b,m,n] = unique(data(nonmissing,[component_mat(i,j) component_mat(i,j+1)]),'rows');
                [b,m,n] = unique(data(:,[component_mat(i,j) component_mat(i,j+1)]),'rows');
                % data_separator(nonmissing,k) = n;
                data_separator(:,k) = n;
                % noalle_separator(k) = size(unique(data(nonmissing,component_mat(i,j))),1)*...
                %    size(unique(data(nonmissing,component_mat(i,j+1))),1);
                noalle_separator(k) = size(unique(data(:,component_mat(i,j))),1)*...
                    size(unique(data(:,component_mat(i,j+1))),1);
                k = k+1;
            end
        end
    end
    
    
elseif strcmp(linkage_model,'independent')
    data_clique = data(:,[1:end-1]);
    data_separator = [];
    noalle_clique = [];
    noalle_separator = [];
    
elseif strcmp(linkage_model,'random')
    adj_mat = mk_adjmat(component_mat);
    nloci = size(adj_mat,1);
    ns = 4*ones(1,nloci);
    porder = [];
    stages = {[1:nloci]};
    clusters = {};
    
    [jtree,root2,cliques,B,w]=graph_to_jtree(adj_mat,ns,porder,stages,clusters);
    C = length(cliques);
    [is,js] = find(jtree > 0);
    separator = cell(C,C);
    for k=1:length(is)
        i = is(k); j = js(k);
        separator{i,j} = find(B(i,:) & B(j,:)); % intersect(cliques{i}, cliques{j});
    end
    
    
    nindividuals = size(data,1);
    if nloci ~= size(data,2)-1 error('Error in transform3.m');
        return
    end
    ncliques = length(cliques);
    [i,j] = find(~cellfun('isempty',triu(separator)));
    ij = [i j];
    nseparators = length(i);
    
    data_clique = zeros(nindividuals, ncliques);
    noalle_clique = zeros(ncliques,1);
    data_separator = zeros(nindividuals, nseparators);
    noalle_separator = zeros(nseparators,1);
    
    k = 1;
    for i = 1:ncliques
        [b,m,n] = unique(data(:,cliques{i}),'rows');
        data_clique(:,k) = n;
        for j = 1:length(cliques{i})
            sz(j) = size(unique(data(:,cliques{i}(j))),1);
        end
        noalle_clique(k) = prod(sz);
        k = k+1;
        clear sz;
    end    
    k = 1;
    for i = 1:nseparators
        [b,m,n] = unique(data(:,separator{ij(i,1),ij(i,2)}),'rows');
        data_separator(:,k) = n;
        
        for j = 1:length(separator{ij(i,1),ij(i,2)})
            sz(j) = size(unique(data(:,separator{ij(i,1),ij(i,2)}(j))),1);
        end
        noalle_separator(k) = prod(sz);
        k = k+1;
        clear sz;
    end
end
clear data;





