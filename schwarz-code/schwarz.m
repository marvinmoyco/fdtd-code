% Schwarz Alternating method 

% Initialization of constants 
l_subdom = 500;               % length of domain
subdom = 4;           % number of subdomains

% initializing of computational domain 
D = randi([0 l_subdom],1,l_subdom);                 % temporary

% initializing of subdomains 
t = round(l_subdom/subdom);

for index = 1:subdom
    x = num2str(index);
    y = strcat('omg',x);

    s.(y) = [];  

    if index == 1
        s.(y) = D(1:t);
    elseif index < subdom
        s.(y) = D(((index-1)*t + 1):index*t);
    else
        s.(y) = D(((index-1)*t + 1):l_subdom);
    end 

end




