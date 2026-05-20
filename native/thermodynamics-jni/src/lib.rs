use jni::objects::JClass;
use jni::sys::{jdouble, jint};
use jni::JNIEnv;

const IDEAL_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

#[no_mangle]
pub extern "system" fn Java_com_example_examplemod_common_rust_ThermodynamicsNative_nativeIdealGasPressure(
    _: JNIEnv,
    _: JClass,
    moles: jdouble,
    temperature_kelvin: jdouble,
    volume_cubic_meters: jdouble,
) -> jdouble {
    if volume_cubic_meters == 0.0 {
        return f64::NAN;
    }

    (moles * IDEAL_GAS_CONSTANT * temperature_kelvin) / volume_cubic_meters
}

#[no_mangle]
pub extern "system" fn Java_com_example_examplemod_common_rust_ThermodynamicsNative_nativeAbiVersion(
    _: JNIEnv,
    _: JClass,
) -> jint {
    1
}

#[cfg(test)]
mod tests {
    use super::IDEAL_GAS_CONSTANT;

    #[test]
    fn ideal_gas_pressure_matches_reference_state() {
        let pressure = (1.0 * IDEAL_GAS_CONSTANT * 273.15) / 0.022_414;
        assert!((pressure - 101_326.0).abs() < 150.0);
    }
}
