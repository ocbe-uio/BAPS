%Seuraavat kolme funktiota liittyvat alkupartition muodostamiseen.

function initial_partition=admixture_initialization(data_matrix,nclusters, Z)

  size_data=size(data_matrix);
  nloci=size_data(2)-1;
  n=max(data_matrix(:,end));
  T=cluster_own(Z,nclusters);
  initial_partition=zeros(size_data(1),1);
  for i=1:n
    kori=T(i);
    here=find(data_matrix(:,end)==i);
    for j=1:length(here)
      initial_partition(here(j),1)=kori;
    end
  end
end
