#!/bin/bash

# Prompt for side to sieve
read -p "Sieve algebraic (a) or rational (r) side? [a/r]: " side_choice
case "$side_choice" in
    r|R)
        side_flag="-r"
        side_name="rational"
        ;;
    a|A|*)
        side_flag="-a"
        side_name="algebraic"
        ;;
esac

echo "Sieving $side_name side"
echo ""

# Prompt for siever version
read -p "Choose siever version [14/15/16]: " siever_choice
case "$siever_choice" in
    14)
        siever="gnfs-lasieve4I14e"
        ;;
    16)
        siever="gnfs-lasieve4I16e"
        ;;
    15|*)
        siever="gnfs-lasieve4I15e"
        ;;
esac

echo "Using $siever"
echo ""

# Prompt for parameters
read -p "Enter start value (in millions): " start_m
read -p "Enter end value (in millions): " end_m

# If start equals end, skip interval prompt
if [ "$start_m" -eq "$end_m" ]; then
    interval_m=1  # Doesn't matter, will only run once
else
    read -p "Enter interval (in millions): " interval_m
fi

# Convert to actual values
start=$((start_m * 1000000))
end=$((end_m * 1000000))
interval=$((interval_m * 1000000))

# Base template file
template="input.job"

# Check if template exists
if [ ! -f "$template" ]; then
    echo "Error: $template not found!"
    exit 1
fi

# Hardcode qintsize
qintsize=1000

# Arrays to store results
declare -a q0_array
declare -a expq_array
declare -a specq_array
declare -a yield_array
declare -a nyld_array
declare -a speed_array

# Iterate through the range
current=$start
index=0
while [ $current -le $end ]; do
    # Create unique filenames
    current_m=$((current / 1000000))
    output_file="output.${current_m}M.out"

    echo ""
    echo "=========================================="
    echo "  Test Sieve Run: ${current_m}M"
    echo "=========================================="

    # Run sieve and capture output
    log_file="log.${current_m}M.txt"
    ./$siever -v -n0 -c $qintsize -f "$current" -o "$output_file" $side_flag "$template" 2>&1 | tee "$log_file"

    # Extract and display summary information
    echo ""
    echo "--- Summary for ${current_m}M ---"

    # Extract values
    specq=$(grep "Special q," "$log_file" | tail -1 | grep -oE "[0-9]+" | head -1)
    yield=$(grep "^Total yield:" "$log_file" | grep -oE "[0-9]+")
    speed=$(grep "sec/rel" "$log_file" | grep -oE "[0-9.]+ sec/rel" | head -1)

    echo "$specq Special q"
    echo "Total yield: $yield"
    echo "($speed)"
    echo ""

    # Calculate expq = qintsize / ln(q0)
    expq=$(echo "scale=3; $qintsize / l($current)" | bc -l)

    # Calculate n-yld = yield * (expq / specq)
    if [ "$specq" -gt 0 ]; then
        nyld=$(echo "scale=6; $yield * ($expq / $specq)" | bc -l)
    else
        nyld="0"
    fi

    # Store in arrays
    q0_array[$index]=$current
    expq_array[$index]=$expq
    specq_array[$index]=$specq
    yield_array[$index]=$yield
    nyld_array[$index]=$nyld
    speed_array[$index]="$speed"

    # Move to next interval
    current=$((current + interval))
    index=$((index + 1))
done

echo ""
echo "=========================================="
echo "All test sieves completed!"
echo "=========================================="
echo ""

# Print formatted table
printf "%-12s %-10s %-8s %-8s %-8s %-8s %-15s %-12s\n" \
    "q0" "qintsize" "expq" "specq" "yield" "n-yld" "speed" "exp_rel"

prev_nyld=""
prev_q0=""

for i in "${!q0_array[@]}"; do
    q0="${q0_array[$i]}"
    expq="${expq_array[$i]}"
    specq="${specq_array[$i]}"
    yield="${yield_array[$i]}"
    nyld="${nyld_array[$i]}"
    speed="${speed_array[$i]}"

    # Calculate exp_rel = avg(last two n-yld) * (q0[i] - q0[i-1]) / qintsize
    if [ -n "$prev_nyld" ]; then
        avg_nyld=$(echo "scale=6; ($nyld + $prev_nyld) / 2" | bc -l)
        q0_diff=$((q0 - prev_q0))
        exp_rel=$(echo "scale=0; $avg_nyld * $q0_diff / $qintsize" | bc -l)
    else
        exp_rel=""
    fi

    # Format values for display
    expq_fmt=$(printf "%.1f" "$expq")
    nyld_fmt=$(printf "%.0f" "$nyld")
    speed_fmt=$(printf "%.3f sec/rel" "${speed% sec/rel}")

    # Print row with proper formatting
    if [ -n "$exp_rel" ]; then
        exp_rel_fmt=$(printf "%.0f" "$exp_rel")
        printf "%-12s %-10s %-8s %-8s %-8s %-8s %-15s %-12s\n" \
            "$q0" "$qintsize" "$expq_fmt" "$specq" "$yield" "$nyld_fmt" "$speed_fmt" "$exp_rel_fmt"
    else
        printf "%-12s %-10s %-8s %-8s %-8s %-8s %-15s\n" \
            "$q0" "$qintsize" "$expq_fmt" "$specq" "$yield" "$nyld_fmt" "$speed_fmt"
    fi

    # Save for next iteration
    prev_nyld="$nyld"
    prev_q0="$q0"
done

echo ""
echo "--- Test Configuration ---"
echo "Siever: $siever"
echo "Side: $side_name"
echo ""

# Calculate suggested lambdas
mfbr=$(grep "^mfbr:" "$template" | awk '{print $2}')
mfba=$(grep "^mfba:" "$template" | awk '{print $2}')
rlim=$(grep "^rlim:" "$template" | awk '{print $2}')
alim=$(grep "^alim:" "$template" | awk '{print $2}')

if [ -n "$mfbr" ] && [ -n "$rlim" ]; then
    rlambda_suggested=$(echo "scale=2; $mfbr / (l($rlim) / l(2))" | bc -l)
    echo "Suggested rlambda: $rlambda_suggested (mfbr=$mfbr / log2(rlim=$rlim))"
fi

if [ -n "$mfba" ] && [ -n "$alim" ]; then
    alambda_suggested=$(echo "scale=2; $mfba / (l($alim) / l(2))" | bc -l)
    echo "Suggested alambda: $alambda_suggested (mfba=$mfba / log2(alim=$alim))"
fi

echo ""
echo "--- Shortened Table ---"
printf "%-12s %-8s %-8s %-12s %-15s\n" "q0" "yield" "n-yld" "exp_rel" "speed"

prev_nyld=""
prev_q0=""
total_exp_rel=0

for i in "${!q0_array[@]}"; do
    q0="${q0_array[$i]}"
    yield="${yield_array[$i]}"
    nyld="${nyld_array[$i]}"
    speed="${speed_array[$i]}"
    expq="${expq_array[$i]}"
    specq="${specq_array[$i]}"

    # Calculate speed in rel/sec = 1/speed
    speed_val="${speed% sec/rel}"
    if [ "$(echo "$speed_val > 0" | bc -l)" -eq 1 ]; then
        speed_relsec=$(echo "scale=3; 1 / $speed_val" | bc -l)
    else
        speed_relsec="0"
    fi

    # Calculate exp_rel = avg(last two n-yld) * (q0[i] - q0[i-1]) / qintsize
    if [ -n "$prev_nyld" ]; then
        avg_nyld=$(echo "scale=6; ($nyld + $prev_nyld) / 2" | bc -l)
        q0_diff=$((q0 - prev_q0))
        exp_rel=$(echo "scale=0; $avg_nyld * $q0_diff / $qintsize" | bc -l)
        total_exp_rel=$(echo "scale=0; $total_exp_rel + $exp_rel" | bc -l)
    else
        exp_rel=""
    fi

    # Format values for display
    nyld_fmt=$(printf "%.0f" "$nyld")
    speed_fmt=$(printf "%.3f rel/sec" "$speed_relsec")

    # Print row with proper formatting
    if [ -n "$exp_rel" ]; then
        exp_rel_fmt=$(printf "%.0f" "$exp_rel")
        printf "%-12s %-8s %-8s %-12s %-15s\n" "$q0" "$yield" "$nyld_fmt" "$exp_rel_fmt" "$speed_fmt"
    else
        printf "%-12s %-8s %-8s %-12s %-15s\n" "$q0" "$yield" "$nyld_fmt" "" "$speed_fmt"
    fi

    # Save for next iteration
    prev_nyld="$nyld"
    prev_q0="$q0"
done

# Print total expected relations
total_exp_rel_fmt=$(printf "%.0f" "$total_exp_rel")
printf "%-12s %-8s %-8s %-12s\n" "" "" "Total:" "$total_exp_rel_fmt"
echo ""

# Create result.job file
cp input.job result.job
if [ "$side_name" = "algebraic" ]; then
    # Ensure the file ends with a newline before appending
    [ -n "$(tail -c1 result.job)" ] && echo "" >> result.job
    echo "lss: 0" >> result.job
fi

# Display result.job contents
echo "--- result.job ---"
cat result.job
echo ""

# Cleanup temporary files
rm -f log.* output.*
