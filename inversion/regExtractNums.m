function nums = regExtractNums(str) 
[a,b] = regexp(str, '\d+(\.\d+)?'); 
nums = zeros(length(a),1); 
for k = 1:length(a) 
    nums(k) = str2double(str(a(k):b(k))); 
end