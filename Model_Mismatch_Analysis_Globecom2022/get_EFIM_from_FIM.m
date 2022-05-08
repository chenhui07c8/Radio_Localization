function EFIM = get_EFIM_from_FIM(FIM, N_states)

    if length(N_states) == 1
        F1 = FIM(1:N_states, 1:N_states);
        F2 = FIM(1:N_states, (N_states+1):end);
        F4 = FIM((N_states+1):end, (N_states+1):end);
        EFIM = F1 - F2*F4^-1*F2.';
    else
        F1 = FIM(1:(N_states(1)-1), 1:(N_states(1)-1));
        F2 = FIM(1:(N_states(1)-1), (N_states(1)):end);
        F4 = FIM((N_states(1)):end, (N_states(1)):end);
%         EFIM = F1 - F2*F4^-1*F2.';
        EFIM = F4 - F2.'*F1^-1*F2;
        if N_states(end)~= size(FIM,1)
            N_states = N_states(end) - N_states(1) + 1;
            F1 = EFIM(1:N_states, 1:N_states);
            F2 = EFIM(1:N_states, (N_states+1):end);
            F4 = EFIM((N_states+1):end, (N_states+1):end);
            EFIM = F1 - F2*F4^-1*F2.';
        end
    end
end