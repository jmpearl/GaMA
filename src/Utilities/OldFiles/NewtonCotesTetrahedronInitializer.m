function [weight] = NewtonCotesTetrahedronInitializer(n)


% symbolic variables 
syms u v w;
phi = {u,v,w,1-u-w-v};

% initialize
i = 0;
j = 0;
k = 0;
l = n;
iter = 1;

while i <= n
    
    q_coord(iter,1:4) = [i,j,k,l]; % coordinate relative to vertices
    
    
    % Calculate Weights
    BF = 1;
    for jj = 1:4
        for ii = 1:q_coord(iter,jj)
            
            BF = BF * (sum(q_coord(iter,:))*phi{jj}-ii+1)/ii;
            
        end
    end


    weight(iter) = double(int(int(int(BF,w,[0,1-u-v]),v,[0,1-u]),u,[0,1]));
    
    % iteration logic
    if l>0
        l=l-1;
        k=k+1;      
    elseif i+j+k==n && i+j<n
        k=0;
        j=j+1;
        l=n-(i+j+k);
    elseif i+j==n && i<n
        i=i+1;
        j=0;
        l=n-(i+j+k);
    else
        break
    end
    iter = iter+1;
    
end

weight=6*weight;

end

