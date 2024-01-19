function out = check_plot(string,plot_list)

out = 0;

n = length(plot_list);
for i = 1:n
    if strcmp(plot_list{i},string)
        out = 1;
    end
end