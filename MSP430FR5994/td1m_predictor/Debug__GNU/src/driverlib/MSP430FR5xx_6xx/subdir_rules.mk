################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Each subdirectory must supply rules for building sources it contributes
src/driverlib/MSP430FR5xx_6xx/%.o: ../src/driverlib/MSP430FR5xx_6xx/%.c $(GEN_OPTS) | $(GEN_FILES) $(GEN_MISC_FILES)
	@echo 'Building file: "$<"'
	@echo 'Invoking: GNU Compiler'
	"/Applications/ti/ccs1250/ccs/tools/compiler/gcc.LTS/bin/msp430-elf-gcc-9.3.1" -c -mmcu=msp430fr5994 -mhwmult=f5series -I"/Applications/ti/ccs1250/ccs/ccs_base/msp430/include_gcc" -I"/Users/giovannicherubini/Desktop/Thesis/Code/t1dm_predictor/MSP430FR5994/td1m_predictor" -I"/Applications/ti/ccs1250/ccs/tools/compiler/gcc.LTS/msp430-elf/include" -O0 -g -gdwarf-3 -gstrict-dwarf -Wall -mlarge -mcode-region=none -mdata-region=lower -MMD -MP -MF"src/driverlib/MSP430FR5xx_6xx/$(basename $(<F)).d_raw" -MT"$(@)"  $(GEN_OPTS__FLAG) -o"$@" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '


