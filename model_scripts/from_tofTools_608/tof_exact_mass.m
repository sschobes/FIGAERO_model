function [mz, namEl, numEl] = tof_exact_mass(compound,nominalmz)

% Get exact mass for compounds, supports nested brackets and isotopes.
% Current limitations: The number following brackets or square brackets
% cannot be larger than 999.
%
% in:
% compound   - cell of compounds
% nominalmz  - 1 returns the nominal mass, 0 doesn't
%
% out:
% mz         - row vector of masses
% namEl      - cell containing the names of the elements of the given compounds
% numEl      - cell containing the number of the elements given in namEl
%
% usage:
% mz = tof_exact_mass('[13C]');
% mz = tof_exact_mass('(H2SO4)30HSO4-', 1);
% [mz, namEl, numEl] = tof_exact_mass('(H2O)6Br-');
%
% Gustaf L?nn
% May 2011

global masstable
elsOut = 0;
electron = 0.0005485799;
isoPresent = 0;

if nargin == 1
    nominalmz = 0;
end

if isempty(masstable)
    load('masstable.mat');
end

masses = masstable.masses;

if ~iscell(compound)
    compound = {compound};
end

if nargout > 1
    elsOut = 1;
    namEl = cell(size(compound));
    numEl = namEl;
end

mz = zeros(size(compound));
oldComp = [];

c = 0;
while c < length(compound)
    c = c + 1;
    comp = compound{c};
    if ~isempty(comp)
        if ~isnumeric(comp)
            if elsOut
                elems = sparse([],[],[],93987,1,0);
            end
            spaces = (comp == ' ');
            comp(spaces) = [];
            
            tempind = 1:length(comp);
            
            % check charge
            if comp(end) == 45
                comp(end) = [];
                if ~nominalmz
                    mz(c) = mz(c) + electron;
                end
            elseif comp(end) == 43
                comp(end) = [];
                if ~nominalmz
                    mz(c) = mz(c) - electron;
                end
            end
            
            % parse brackets
            leftb = tempind(comp == 40);
            rightb = tempind(comp == 41);
            if length(leftb) ~= 0
                i = 0;
                while i < length(leftb)
                    i = i + 1;
                    weight = 1;
                    % determine the correct index for rightb if nested
                    % brackets present
                    rbInd = i;
                    inc2 = 0;
                    if i ~= length(leftb)
                       inc1 = 1;
                       while leftb(i + inc1) < rightb(i + inc2)
                           inc1 = inc1 + 1;
                           inc2 = inc2 + 1;
                           if i + inc1 > length(leftb)
                               rbInd = i + inc2;
                               break;
                           end
                           rbInd = i + inc2;
                       end
                    end
                    m = min(3,length(comp)-rightb(rbInd));
                    p = comp(rightb(rbInd)+1:rightb(rbInd)+m);
                    nums = p > 47 & p < 59;
                    first0 = find(nums == 0);
                    if length(first0) ~= 0
                        first0 = first0(1);
                        nums = nums(1:first0-1);
                    end
                    t = tempind(nums);
                    if isempty(t) || nums(1) == 0
                        t = 0;
                        nums = [];
                    end
                    if any(nums)
                        weight = Gstr2int(p(nums));
                    end
                    if elsOut
                        newComp = comp(leftb(i)+1:rightb(rbInd)-1);
                        varNewComp = strrep(newComp, '(', 'j'); % turn into valid field name, jkpq not used in compounds
                        varNewComp = strrep(varNewComp, ')', 'k');
                        varNewComp = strrep(varNewComp, '[', 'p');
                        varNewComp = strrep(varNewComp, ']', 'q');
                        if isfield(oldComp, varNewComp)
                            mzT = oldComp.(varNewComp).mz;
                            namT = oldComp.(varNewComp).nam;
                            numT = oldComp.(varNewComp).num;
                        else
                            [mzT, namT, numT] = tof_exact_mass(newComp,nominalmz);
                            try
                                oldComp.(varNewComp).mz = mzT;
                                oldComp.(varNewComp).nam = namT;
                                oldComp.(varNewComp).num = numT;
                            end
                        end
                        mz(c) = mz(c) + mzT * weight;
                        for j = 1:length(namT{1})
                            iso = 0;
                            lnamT = length(namT{1}{j});
                            if lnamT > 2
                                numsl = namT{1}{j} > 47 & namT{1}{j} < 59;
                                numsl = namT{1}{j}(numsl);
                                iso = Gstr2int(numsl);
                                namT{1}{j} = namT{1}{j}(length(numsl)+2:end-1);
                                lnamT = length(namT{1}{j});
                            end
                            if lnamT > 1
                                ind = (namT{1}{j}(1) - 64) * 59^2 + (namT{1}{j}(2) - 64) * 59 + iso;
                            else
                                ind = (namT{1}{j} - 64) * 59^2 + iso;
                            end
                            elems(ind) = elems(ind) + numT{1}(j) * weight;
                        end
                    else
                        newComp = comp(leftb(i)+1:rightb(rbInd)-1);
                        varNewComp = strrep(newComp, '(', 'j'); % turn into valid field name, jkpq not used in compounds
                        varNewComp = strrep(varNewComp, ')', 'k');
                        varNewComp = strrep(varNewComp, '[', 'p');
                        varNewComp = strrep(varNewComp, ']', 'q');
                        if isfield(oldComp, varNewComp)
                            mz(c) = mz(c) + oldComp.(varNewComp) * weight;
                        else
                            newMass = tof_exact_mass(newComp,nominalmz);
                            try
                                oldComp.(varNewComp) = newMass;
                            end
                            mz(c) = mz(c) + newMass * weight;
                        end
                    end
                    comp(leftb(i):rightb(rbInd)+t(end)) = ' ';
                    if inc2 ~= 0
                        i = i + inc2;
                    end
                end
                comp(comp == 32) = [];
            end
            
            % parse square brackets & calculate isotope masses
            leftb = tempind(comp == 91);
            rightb = tempind(comp == 93);
            if length(leftb) ~= 0
                isoPresent = 1;
                for i = 1:length(leftb)
                    weight = 1;
                    p = comp(leftb(i)+1:rightb(i)-1);
                    nums = p > 47 & p < 59;
                    iso = Gstr2int(p(nums));
                    p = p(~nums);
                    ind = fix(p);
                    if length(ind) > 1
                        ind = (ind(1) - 64) * 59^2 + (ind(2) - 64) * 59 + iso;
                    else
                        ind = (ind - 64) * 59^2 + iso;
                    end
                    m = min(3,length(comp)-rightb(i));
                    pr = comp(rightb(i)+1:rightb(i)+m);
                    nums = pr > 47 & pr < 59;
                    first0 = find(nums == 0);
                    if length(first0) ~= 0
                        first0 = first0(1);
                        nums = nums(1:first0-1);
                    end
                    t = tempind(nums);
                    if isempty(t) || nums(1) == 0
                        t = 0;
                        nums = [];
                    end
                    if any(nums)
                        weight = Gstr2int(pr(nums));
                    end
                    if elsOut
                        elems(ind) = elems(ind) + weight;
                    end
                    addm = masses(ind);
                    if addm ~= 0
                        if nominalmz
                            mz(c) = mz(c) + round(addm) * weight;
                        else
                            mz(c) = mz(c) + addm * weight;
                        end
                    else
                        disp(['tof_exact_mass: Isotope ' num2str(iso) p ' not found in database.']);
                        mz(c) = NaN;
                    end
                    comp(leftb(i):rightb(i)+t(end)) = ' ';
                end
                comp(comp == 32) = [];
            end
            
            % calculate standard masses
            parts = tempind(comp > 64 & comp < 91);
            parts(end + 1) = length(comp) + 1;
            
            for i = 1:length(parts) - 1
                weight = 1;
                p = comp(parts(i):parts(i+1)-1);
                nums = p > 47 & p < 59;
                if any(nums)
                    weight = Gstr2int(p(nums));
                    p = p(~nums);
                end
                ind = fix(p);
                if length(ind) > 1
                    ind = (ind(1) - 64) * 59^2 + (ind(2) - 64) * 59;
                else
                    ind = (ind - 64) * 59^2;
                end
                addm = masses(ind);
                if elsOut
                    elems(ind) = elems(ind) + weight;
                end
                if addm ~= 0
                    if nominalmz
                        mz(c) = mz(c) + round(addm) * weight;
                    else
                        mz(c) = mz(c) + addm * weight;
                    end
                else
                    disp(['tof_exact_mass: Compound ' p ' not found in database.']);
                    mz(c) = NaN;
                end
            end
            
            % process elements
            if elsOut && ~isnan(mz(c))
                inds = find(elems ~= 0);
                linds = length(inds);
                namEl{c} = cell(1,linds);
                numEl{c} = zeros(1,linds);
                for i = 1:linds
                    el = masstable.names{masstable.index(inds(i))};
                    namEl{c}{i} = el;
                    numEl{c}(i) = elems(inds(i));
                end
                if isoPresent
                    elems = sparse([],[],[],93987,1,0);
                    for j = 1:length(namEl{c})
                        iso = 0;
                        lnamT = length(namEl{c}{j});
                        if lnamT > 2
                            numsl = namEl{c}{j} > 47 & namEl{c}{j} < 59;
                            numsl = namEl{c}{j}(numsl);
                            iso = Gstr2int(numsl);
                            namEl{c}{j} = namEl{c}{j}(length(numsl)+2:end-1);
                            lnamT = length(namEl{c}{j});
                        end
                        if lnamT > 1
                            ind = (namEl{c}{j}(1) - 64) * 59^2 + (namEl{c}{j}(2) - 64) * 59 + iso;
                        else
                            ind = (namEl{c}{j} - 64) * 59^2 + iso;
                        end
                        elems(ind) = elems(ind) + numEl{c}(j);
                    end
                    inds = find(elems ~= 0);
                    linds = length(inds);
                    namEl{c} = cell(1,linds);
                    numEl{c} = zeros(1,linds);
                    for i = 1:linds
                        el = masstable.names{masstable.index(inds(i))};
                        namEl{c}{i} = el;
                        numEl{c}(i) = elems(inds(i));
                    end
                end
            end
        else
            if length(comp) > 0
                mz(c) = comp;
            else
                mz(c) = NaN;
            end
            if elsOut
                namEl{c} = [];
                numEl{c} = [];
            end
        end
    else
        mz(c)=NaN;
        if elsOut
            namEl{c} = [];
            numEl{c} = [];
        end
    end
end


function num = Gstr2int(str)

num = 0;
t = length(str);
for i = 1:t
    num = num + (str(i) - 48) * 10^(t-i);
end