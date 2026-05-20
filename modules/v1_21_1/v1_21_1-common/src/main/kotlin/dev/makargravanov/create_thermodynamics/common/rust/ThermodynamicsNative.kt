package dev.makargravanov.create_thermodynamics.common.rust

object ThermodynamicsNative {
    init {
        NativeLibraryLoader.load()
    }

    fun idealGasPressure(moles: Double, temperatureKelvin: Double, volumeCubicMeters: Double): Double =
        nativeIdealGasPressure(moles, temperatureKelvin, volumeCubicMeters)

    fun abiVersion(): Int = nativeAbiVersion()

    @JvmStatic
    private external fun nativeIdealGasPressure(
        moles: Double,
        temperatureKelvin: Double,
        volumeCubicMeters: Double,
    ): Double

    @JvmStatic
    private external fun nativeAbiVersion(): Int
}
