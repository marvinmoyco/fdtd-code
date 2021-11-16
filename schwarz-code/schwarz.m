% Schwarz Alternating method 

% Initialization of constants 
l_dom = 500;                                 % length of domain      % arbitrary 
subdom = 4;                                  % number of subdomains  % arbitrary 

% initializing of computational domain 
D = randi([0 l_dom],1,l_dom);                % temporary

% initializing of subdomains 
t = round(l_dom/subdom);

for index = 1:subdom
    x = num2str(index);
    y = strcat('omg',x);

    omegas.(y) = [];  

    if index == 1
        omegas.(y) = D(1:t);
    elseif index < subdom
        omegas.(y) = D(((index-1)*t + 1):index*t);
    else
        omegas.(y) = D(((index-1)*t + 1):l_dom);
    end 

end

% set external boundary conditions 

% solve PDE on subdomain
% insert 1D FDTD 

% placeholder
for index = 1:subdom
    x = num2str(index);
    y = strcat('omg',x);
    z = strcat('omg_slvd',x);

    omegas_slvd.(z) = omegas.(y) + 1;
end  

% Transfer of data 
olap_sz = 20; 









