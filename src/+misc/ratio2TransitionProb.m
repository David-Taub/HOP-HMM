function P = ratio2TransitionProb(lengthOfB, AtoBRatio)
    P = AtoBRatio / lengthOfB * (1 - AtoBRatio);
end