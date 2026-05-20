package com.example.examplemod.common

import com.example.examplemod.common.rust.ThermodynamicsNative

object CommonModule {
    fun demoPressureKilopascals(): Double =
        ThermodynamicsNative.idealGasPressure(
            moles = 1.0,
            temperatureKelvin = 273.15,
            volumeCubicMeters = 0.022414,
        ) / 1_000.0

    fun nativeAbiVersion(): Int = ThermodynamicsNative.abiVersion()
}
