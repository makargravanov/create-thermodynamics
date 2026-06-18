use jni::objects::{JClass, JDoubleArray, JObject, JObjectArray, JString};
use jni::sys::{jboolean, jdouble, jint, jlong, jobjectArray, JNI_FALSE, JNI_TRUE};
use jni::JNIEnv;

pub mod chemistry;

const IDEAL_GAS_CONSTANT: f64 = 8.314_462_618_153_24;

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeIdealGasPressure(
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
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeAbiVersion(
    _: JNIEnv,
    _: JClass,
) -> jint {
    4
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeReplaceMinecraftItemChemicalBindings(
    mut env: JNIEnv,
    _: JClass,
    item_ids: JObjectArray,
    substance_ids: JObjectArray,
    mol_per_items: JDoubleArray,
) {
    match read_item_chemical_bindings_from_jvm(&mut env, item_ids, substance_ids, mol_per_items)
        .and_then(chemistry::minecraft::chem_api::replace_item_chemical_bindings)
    {
        Ok(()) => {}
        Err(error) => throw_java_exception(
            &mut env,
            "java/lang/IllegalArgumentException",
            &error.to_string(),
        ),
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeClearMinecraftItemChemicalBindings(
    mut env: JNIEnv,
    _: JClass,
) {
    if let Err(error) = chemistry::minecraft::chem_api::clear_item_chemical_bindings() {
        throw_java_exception(
            &mut env,
            "java/lang/IllegalStateException",
            &error.to_string(),
        );
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeMinecraftItemChemicalBindingCount(
    mut env: JNIEnv,
    _: JClass,
) -> jint {
    match chemistry::minecraft::chem_api::item_chemical_binding_count() {
        Ok(count) if count <= jint::MAX as usize => count as jint,
        Ok(count) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalStateException",
                &format!("minecraft item chemical binding count {count} does not fit into Int"),
            );
            -1
        }
        Err(error) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalStateException",
                &error.to_string(),
            );
            -1
        }
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeHasMinecraftItemChemicalBinding(
    mut env: JNIEnv,
    _: JClass,
    item_id: JString,
) -> jboolean {
    let result = read_java_string(&mut env, item_id)
        .and_then(|item_id| chemistry::minecraft::chem_api::has_item_chemical_binding(&item_id));
    match result {
        Ok(true) => JNI_TRUE,
        Ok(false) => JNI_FALSE,
        Err(error) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalStateException",
                &error.to_string(),
            );
            JNI_FALSE
        }
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeStaticSubstanceIds(
    mut env: JNIEnv,
    _: JClass,
) -> jobjectArray {
    match chemistry::minecraft::chem_api::static_substance_ids()
        .and_then(|ids| java_string_array(&mut env, &ids))
    {
        Ok(array) => array,
        Err(error) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalStateException",
                &error.to_string(),
            );
            std::ptr::null_mut()
        }
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeCreateSingleZoneReactor(
    mut env: JNIEnv,
    _: JClass,
    volume_cubic_meters: jdouble,
    item_input_count: jint,
    item_output_count: jint,
    fluid_input_count: jint,
    fluid_output_count: jint,
) -> jlong {
    let result = read_non_negative_count("itemInputCount", item_input_count)
        .and_then(|item_inputs| {
            Ok((
                item_inputs,
                read_non_negative_count("itemOutputCount", item_output_count)?,
                read_non_negative_count("fluidInputCount", fluid_input_count)?,
                read_non_negative_count("fluidOutputCount", fluid_output_count)?,
            ))
        })
        .and_then(|(item_inputs, item_outputs, fluid_inputs, fluid_outputs)| {
            chemistry::minecraft::chem_api::register_single_zone_reactor(
                volume_cubic_meters,
                item_inputs,
                item_outputs,
                fluid_inputs,
                fluid_outputs,
            )
        });
    match result {
        Ok(reactor_id) if reactor_id.0 <= jlong::MAX as u64 => reactor_id.0 as jlong,
        Ok(reactor_id) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalStateException",
                &format!("reactor id {} does not fit into Long", reactor_id.0),
            );
            -1
        }
        Err(error) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalArgumentException",
                &error.to_string(),
            );
            -1
        }
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeRemoveReactor(
    mut env: JNIEnv,
    _: JClass,
    reactor_id: jlong,
) {
    let result =
        read_reactor_id(reactor_id).and_then(chemistry::minecraft::chem_api::remove_reactor);
    if let Err(error) = result {
        throw_java_exception(
            &mut env,
            "java/lang/IllegalArgumentException",
            &error.to_string(),
        );
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeReactorCount(
    mut env: JNIEnv,
    _: JClass,
) -> jint {
    match chemistry::minecraft::chem_api::reactor_count() {
        Ok(count) if count <= jint::MAX as usize => count as jint,
        Ok(count) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalStateException",
                &format!("reactor count {count} does not fit into Int"),
            );
            -1
        }
        Err(error) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalStateException",
                &error.to_string(),
            );
            -1
        }
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeTickReactor(
    mut env: JNIEnv,
    _: JClass,
    reactor_id: jlong,
    dt_seconds: jdouble,
) {
    let result = read_reactor_id(reactor_id).and_then(|reactor_id| {
        chemistry::minecraft::chem_api::tick_reactor(reactor_id, dt_seconds)
    });
    if let Err(error) = result {
        throw_java_exception(
            &mut env,
            "java/lang/IllegalArgumentException",
            &error.to_string(),
        );
    }
}

#[no_mangle]
pub extern "system" fn Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeInsertItemStackToReactorInput(
    mut env: JNIEnv,
    _: JClass,
    reactor_id: jlong,
    input_index: jint,
    item_id: JString,
    item_count: jint,
) -> jdouble {
    let result = read_reactor_id(reactor_id)
        .and_then(|reactor_id| {
            Ok((
                reactor_id,
                read_non_negative_count("inputIndex", input_index)?,
                read_positive_u32("itemCount", item_count)?,
                read_java_string(&mut env, item_id)?,
            ))
        })
        .and_then(|(reactor_id, input_index, item_count, item_id)| {
            chemistry::minecraft::chem_api::insert_item_stack_to_reactor_input(
                reactor_id,
                input_index,
                &item_id,
                item_count,
            )
        });
    match result {
        Ok(mol_inserted) => mol_inserted,
        Err(error) => {
            throw_java_exception(
                &mut env,
                "java/lang/IllegalArgumentException",
                &error.to_string(),
            );
            f64::NAN
        }
    }
}

fn read_item_chemical_bindings_from_jvm(
    env: &mut JNIEnv,
    item_ids: JObjectArray,
    substance_ids: JObjectArray,
    mol_per_items: JDoubleArray,
) -> chemistry::ChemistryResult<Vec<chemistry::minecraft::chem_api::ItemChemicalBinding>> {
    let item_count = env
        .get_array_length(&item_ids)
        .map_err(|error| jni_error_to_chemistry_error("itemIds", error))?;
    let substance_count = env
        .get_array_length(&substance_ids)
        .map_err(|error| jni_error_to_chemistry_error("substanceIds", error))?;
    let amount_count = env
        .get_array_length(&mol_per_items)
        .map_err(|error| jni_error_to_chemistry_error("molPerItems", error))?;
    if item_count != substance_count || item_count != amount_count {
        return Err(chemistry::ChemistryError::InvalidMixtureState(format!(
            "itemIds, substanceIds and molPerItems lengths must match, got {item_count}, {substance_count}, {amount_count}"
        )));
    }

    let mut amounts = vec![0.0; amount_count as usize];
    env.get_double_array_region(&mol_per_items, 0, &mut amounts)
        .map_err(|error| jni_error_to_chemistry_error("molPerItems", error))?;

    let mut bindings = Vec::with_capacity(item_count as usize);
    for index in 0..item_count {
        let item_id = {
            let value = env
                .get_object_array_element(&item_ids, index)
                .map_err(|error| jni_error_to_chemistry_error("itemIds", error))?;
            read_java_string(env, JString::from(value))?
        };
        let substance_id = {
            let value = env
                .get_object_array_element(&substance_ids, index)
                .map_err(|error| jni_error_to_chemistry_error("substanceIds", error))?;
            read_java_string(env, JString::from(value))?
        };
        bindings.push(chemistry::minecraft::chem_api::ItemChemicalBinding::new(
            item_id,
            substance_id.as_str(),
            amounts[index as usize],
        ));
    }

    Ok(bindings)
}

fn read_java_string(env: &mut JNIEnv, value: JString) -> chemistry::ChemistryResult<String> {
    env.get_string(&value)
        .map(|value| value.into())
        .map_err(|error| jni_error_to_chemistry_error("String", error))
}

fn read_non_negative_count(name: &str, value: jint) -> chemistry::ChemistryResult<usize> {
    if value < 0 {
        return Err(chemistry::ChemistryError::InvalidMixtureState(format!(
            "{name} must be non-negative, got {value}"
        )));
    }
    Ok(value as usize)
}

fn read_positive_u32(name: &str, value: jint) -> chemistry::ChemistryResult<u32> {
    if value <= 0 {
        return Err(chemistry::ChemistryError::InvalidMixtureState(format!(
            "{name} must be positive, got {value}"
        )));
    }
    Ok(value as u32)
}

fn read_reactor_id(
    reactor_id: jlong,
) -> chemistry::ChemistryResult<chemistry::minecraft::worker::reactors_worker::ReactorInstanceId> {
    if reactor_id < 0 {
        return Err(chemistry::ChemistryError::InvalidMixtureState(format!(
            "reactorId must be non-negative, got {reactor_id}"
        )));
    }
    Ok(chemistry::minecraft::worker::reactors_worker::ReactorInstanceId(reactor_id as u64))
}

fn jni_error_to_chemistry_error(
    context: &str,
    error: jni::errors::Error,
) -> chemistry::ChemistryError {
    chemistry::ChemistryError::InvalidMixtureState(format!(
        "failed to read {context} from JVM: {error}"
    ))
}

fn throw_java_exception(env: &mut JNIEnv, class_name: &str, message: &str) {
    let _ = env.throw_new(class_name, message);
}

fn java_string_array(
    env: &mut JNIEnv,
    values: &[String],
) -> chemistry::ChemistryResult<jobjectArray> {
    let string_class = env
        .find_class("java/lang/String")
        .map_err(|error| jni_error_to_chemistry_error("java/lang/String class", error))?;
    let array = env
        .new_object_array(values.len() as jint, string_class, JObject::null())
        .map_err(|error| jni_error_to_chemistry_error("String array", error))?;
    for (index, value) in values.iter().enumerate() {
        let java_value = env
            .new_string(value)
            .map_err(|error| jni_error_to_chemistry_error("String", error))?;
        env.set_object_array_element(&array, index as jint, java_value)
            .map_err(|error| jni_error_to_chemistry_error("String array", error))?;
    }
    Ok(array.into_raw())
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
