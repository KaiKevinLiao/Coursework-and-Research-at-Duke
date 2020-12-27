cd 'C:\Users\kaike\OneDrive\Desktop\Non_Market_Valuation\Problem Set -- Hedonics';

A = readmatrix('Non_Paramatric_Hedonic_Price_Gradient.csv');
plot(A(:,1),A(:,3))
hold on
plot(A(:,1),A(:,5))
plot(A(:,1),A(:,7))
plot(A(:,1),A(:,9))
hold off
legend('h=1','h=3','h=10','h=1000')
title('Task 4 - Hedonic Price Gradient for Different Values of Bandwidth')
xlabel('\chi') 
ylabel('\beta (\chi)')