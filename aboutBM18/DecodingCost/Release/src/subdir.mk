################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BothMay.cpp \
../src/BothMayOld.cpp \
../src/Decoder.cpp \
../src/DecodingCost.cpp \
../src/MayOzerov.cpp \
../src/Prange.cpp \
../src/SternMayOzerov.cpp \
../src/Tools.cpp 

CPP_DEPS += \
./src/BothMay.d \
./src/BothMayOld.d \
./src/Decoder.d \
./src/DecodingCost.d \
./src/MayOzerov.d \
./src/Prange.d \
./src/SternMayOzerov.d \
./src/Tools.d 

OBJS += \
./src/BothMay.o \
./src/BothMayOld.o \
./src/Decoder.o \
./src/DecodingCost.o \
./src/MayOzerov.o \
./src/Prange.o \
./src/SternMayOzerov.o \
./src/Tools.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++2a -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/BothMay.d ./src/BothMay.o ./src/BothMayOld.d ./src/BothMayOld.o ./src/Decoder.d ./src/Decoder.o ./src/DecodingCost.d ./src/DecodingCost.o ./src/MayOzerov.d ./src/MayOzerov.o ./src/Prange.d ./src/Prange.o ./src/SternMayOzerov.d ./src/SternMayOzerov.o ./src/Tools.d ./src/Tools.o

.PHONY: clean-src

