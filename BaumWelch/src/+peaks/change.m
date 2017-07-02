% a = dir;
% b = {a.name};
% b = {b{3:end}};
% oldB = b;

% for i =1:length(b); b{i} = strrep(b{i}, 'H3K27ac-', ''); end;
% for i =1:length(b); b{i} = strrep(b{i}, 'H3K4me1-', ''); end;
% for i =1:length(b); b{i} = strrep(b{i}, 'H3K4me3-', ''); end;
% for i =1:length(b); b{i} = strrep(b{i}, '.peaks', ''); end;
% for i =1:length(b); b{i} = strrep(b{i}, '.txt', ''); end;
% for i =1:length(b); b{i} = strrep(b{i}, 'RenLab', ''); end;
% for i =1:length(b); movefile(oldB{i}, b{i}); end;

function change
    b = dir('/cs/cbio/david/data/mat/H3K27ac');
    b = {b.name};
    b = {b{3:end}};
    for i = 1:length(b)
        b{i}
        load(['/cs/cbio/david/data/mat/H3K27ac/', b{i}]);
        H3K27ac = [S{:}];
        load(['/cs/cbio/david/data/mat/H3K4me1/', b{i}]);
        H3K4me1 = [S{:}];
        load(['/cs/cbio/david/data/mat/H3K4me3/', b{i}]);
        H3K4me3 = [S{:}];
        T(i).name = b{i};
        T(i).H3K27ac.chr        = {H3K27ac.chr};
        T(i).H3K27ac.from       = [H3K27ac.from];
        T(i).H3K27ac.to         = [H3K27ac.to];
        T(i).H3K27ac.height     = [H3K27ac.height];
        T(i).H3K27ac.gTSSdist   = [H3K27ac.gTSSdist];
        T(i).H3K4me1.chr        = {H3K4me1.chr};
        T(i).H3K4me1.from       = [H3K4me1.from];
        T(i).H3K4me1.to         = [H3K4me1.to];
        T(i).H3K4me1.height     = [H3K4me1.height];
        T(i).H3K4me3.chr        = {H3K4me3.chr};
        T(i).H3K4me3.from       = [H3K4me3.from];
        T(i).H3K4me3.to         = [H3K4me3.to];
        T(i).H3K4me3.height     = [H3K4me3.height];
    end
    save('/cs/cbio/david/data/mat/peaks_raw.mat', 'T');
end
