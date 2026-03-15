function S = read_sph(fn)

if nargin < 1
    fn = 'test.sph';
end

#L = strsplit(strrep(fileread(fn), sprintf('\r'), ''), sprintf('\n'));
L = strsplit(strrep(fileread(fn), sprintf('\r'), ''), sprintf('\n'), 'collapsedelimiters', false);

S = struct([]);
k = 0;
i = 1;

while i <= numel(L)-8

    % Header check
    a = sscanf(L{i+2},'%f');
    f = regexp(L{i+3},'Frequency\s*=\s*([0-9eE\+\-\.]+)','tokens','once');

    if numel(a) >= 4 && ~isempty(f)

        k = k + 1;

        ntheta = a(1);
        nphi   = a(2);
        nmax   = a(3);
        mmax   = a(4);
        freq   = str2double(f{1});

        mp = zeros(mmax+1,2);
        A  = complex(zeros(2,2*mmax+1,nmax));

        row = i + 8;

        m = -1;
        n = 0;

        while row <= numel(L)

            nums = sscanf(L{row},'%f');

            if numel(nums) == 2

                m = nums(1);
                mp(m+1,:) = nums';

                if m == 0
                    n = 1;
                else
                    n = m;
                end

            elseif numel(nums) == 4

                if m == 0

                    A(1,mmax+1,n) = complex(nums(1),nums(2));
                    A(2,mmax+1,n) = complex(nums(3),nums(4));
                    n = n + 1;

                else

                    A(1,-m+mmax+1,n) = complex(nums(1),nums(2));
                    A(2,-m+mmax+1,n) = complex(nums(3),nums(4));

                    row = row + 1;
                    nums = sscanf(L{row},'%f');

                    A(1,m+mmax+1,n) = complex(nums(1),nums(2));
                    A(2,m+mmax+1,n) = complex(nums(3),nums(4));

                    n = n + 1;

                end

            else
                break
            end

            row = row + 1;

        end

        S(k).freq   = freq;
        S(k).ntheta = ntheta;
        S(k).nphi   = nphi;
        S(k).nmax   = nmax;
        S(k).mmax   = mmax;
        S(k).mp     = mp;
        S(k).Q_ha   = conj(A)*sqrt(8*pi);
        S(k).Q_ti   = A;
        S(k).P_ha   = 0.5*sum( ( (conj(A)*sqrt(8*pi)) .* conj( conj(A)*sqrt(8*pi) ) )(:) )
        S(k).P_ti   = sum(mp(:,2))*8*pi
        i = row;

    else
        i = i + 1;
    end

end


fprintf('\n');
fprintf('Qmns Hansen');
cellfun(@(x) disp(permute(x,[2 3 1])),{S.Q_ha})
fprintf('\n');
fprintf('Qmns Ticra');
cellfun(@(x) disp(permute(x,[2 3 1])),{S.Q_ti})

fprintf('\n\n\n\n\n\n\n Read Q Files \n ------------\n');

fprintf('\n%3s %8s %8s %8s %8s %14s %18s %18s\n','Nr','nTheta','nPhi','nMax','mMax','Freq in GHz','Total Power Hansen','Total Power Ticra');

for j = 1:numel(S)
    fprintf('%3d %8d %8d %8d %8d %14.6g %18.6g %18.6g\n',...
        j,S(j).ntheta,S(j).nphi,S(j).nmax,S(j).mmax,S(j).freq/1e9,S(j).P_ha,S(j).P_ti);
end
fprintf('\n\n\n\n\n\n\n');
end
