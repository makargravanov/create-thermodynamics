package dev.makargravanov.create_thermodynamics.common

import dev.makargravanov.create_thermodynamics.common.chemistry.binding.DefaultItemChemicalBindings
import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative

object CommonModule {
    fun initializeNativeChemistry() {
        ThermodynamicsNative.configureItemChemicalBindings(DefaultItemChemicalBindings.bindings)
    }

    fun demoPressureKilopascals(): Double =
        ThermodynamicsNative.idealGasPressure(
            moles = 1.0,
            temperatureKelvin = 273.15,
            volumeCubicMeters = 0.022414,
        ) / 1_000.0

    fun nativeAbiVersion(): Int = ThermodynamicsNative.abiVersion()

    fun nativeStaticSubstanceIds(): List<String> = ThermodynamicsNative.staticSubstanceIds()
}
