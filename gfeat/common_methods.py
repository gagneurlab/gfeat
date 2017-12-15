def boyer_moore_search(text, pattern):
    occurrences = dict()
    for letter in {'A', 'C', 'G', 'T'}:
        occurrences[letter] = pattern.rfind(letter)
    m = len(pattern)
    n = len(text)
    i = m - 1  # text index
    j = m - 1  # pattern index
    count = 0;
    while i < n:
        if text[i] == pattern[j]:
            if j == 0:
                count += 1
                i += 2*m - 1
                j = m - 1
            else:
                i -= 1
                j -= 1
        else:
            l = occurrences[text[i]]
            i = i + m - min(j, 1 + l)
            j = m - 1
    return count
