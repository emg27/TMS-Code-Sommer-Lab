%%Clustering plot fixer
function[colorvalue]=plotclusters(linez,currentcluster)

%currentcluster=matrixpt-lowcluster;

switch currentcluster
    case 0
        colorvalue='k';
    case 1
        colorvalue='b';
    case 2
        colorvalue='r';
    case 3
        colorvalue='g';
    case 4
        colorvalue='m';
    case 5
        colorvalue='c';
    case 6
        colorvalue=[.6 .9 .7];
    case 7
        colorvalue=[.5 0 .5];
    case 8
        colorvalue=[.4 .75 0];
    case 9
        colorvalue=[0 .5 .5];
end
set(linez,'Color',colorvalue);
set(linez,'LineWidth',1);
end