import pytest
from src.general.configfile_tools import check_file_path
from src.general.configfile_tools import parse_dir_path
# from src.general.configfile_tools import strip_of_comments
# from src.general.configfile_tools import clean_lines
from src.general.configfile_tools import get_configurations
from src.general.configfile_tools import configurations_consistency_check


def test_check_file_path_file():
    file_path = "tests/test_general/test_configfile_tools.py"
    assert check_file_path(file_path) == file_path


def test_check_file_path_file():
    file_path = "tests/test_general/NO_SUCH_FILE"
    with pytest.raises(FileNotFoundError) as exc_info:
        check_file_path(file_path)
        assert "No such file or directory: " in str(exc_info)


def test_check_file_path_dir():
    file_path = "tests/test_general/"
    with pytest.raises(FileNotFoundError) as exc_info:
        check_file_path(file_path)
        assert "No such file or directory: " in str(exc_info)


def test_parse_dir_path_file():
    dir_path = "tests/test_general/test_configfile_tools.py"
    answer = "tests/test_general/"
    assert parse_dir_path(dir_path) == answer


def test_parse_dir_path_non_existing_dir():
    dir_path = "tests/test_general/NO_SUCH_DIR"
    with pytest.raises(FileNotFoundError) as exc_info:
        parse_dir_path(dir_path)
        assert "No such file or directory: " in str(exc_info)


def test_parse_dir_path_empty_dir():
    dir_path = ""
    with pytest.raises(FileNotFoundError) as exc_info:
        parse_dir_path(dir_path)
        assert "No such file or directory: " in str(exc_info)


def test_parse_dir_current_dir():
    dir_path = "."
    answer = "./"
    assert parse_dir_path(dir_path) == answer


def test_parse_dir_upper_dir():
    dir_path = ".."
    answer = "../"
    assert parse_dir_path(dir_path) == answer


def test_parse_dir_from_upper_dir():
    dir_path = "src/../tests/"
    answer = "src/../tests/"
    print(parse_dir_path(dir_path))
    assert parse_dir_path(dir_path) == answer





def test_get_configurations_check_names_Gen_Mem_Inter_correct():
    example = [
        "-GeneralParameters",
        "ProjectName SM_Ah_GKB_L",
        "SimulationTimeInPs 100000",
        "StepSizeInFs 200",
        "ReportIntervalInFs 100000",
        "Integrator LM",
        "FrictionIninvertPs 0.01",
        "setVelocitiesToTemperature true",
        "EpsilonH 4",
        "TemperatureinKelvin 300",
        "TextOutputs  XYZ PSF",
        "BinOutputs   XYZ TPK",
        "-Membrane  0",
        "MeshFile USphere_1002d_r_s0.ply",
        "NodeMass 50",
        "SpringModel N",
        "MeanCurvatureModel ItzyksonBarycentric",
        "MeanCurvatureCoeff 50",
        "Scale 1000",
        "SurfaceConstraint GWCAH 20",
        "SurfaceConstraintValue Au",
        "SurfaceConstraintCoeff 1.3",
        "VolumeConstraint N",
        "-InteractionTable",
        "M0 0",
    ]
    configs = get_configurations(example)
    assert configs[0].get_name() == "GeneralParameters"
    assert configs[1].get_name() == "Membrane"
    assert configs[2].get_name() == "InteractionTable"


def test_get_configurations_check_names_Gen_Mem_Inter_incorrect():
    example = [
        "-GeneralParametersw",
        "ProjectName SM_Ah_GKB_L",
        "SimulationTimeInPs 100000",
        "StepSizeInFs 200",
        "ReportIntervalInFs 100000",
        "Integrator LM",
        "FrictionIninvertPs 0.01",
        "setVelocitiesToTemperature true",
        "EpsilonH 4",
        "TemperatureinKelvin 300",
        "TextOutputs  XYZ PSF",
        "BinOutputs   XYZ TPK",
        "-Membrane  0",
        "MeshFile USphere_1002d_r_s0.ply",
        "NodeMass 50",
        "SpringModel N",
        "MeanCurvatureModel ItzyksonBarycentric",
        "MeanCurvatureCoeff 50",
        "Scale 1000",
        "SurfaceConstraint GWCAH 20",
        "SurfaceConstraintValue Au",
        "SurfaceConstraintCoeff 1.3",
        "VolumeConstraint N",
        "-InteractionTable",
        "M0 0",
    ]
    with pytest.raises(TypeError) as exc_info:
        configs = get_configurations(example)
        assert (
            ' is not a section name. First non-comment statement of a configuration file must begin with a section declaration.\nSection declarations must be at the beginning of the line and in the following format "-SectionName". SectionName choices are: '
            in str(exc_info)
        )


def test_configurations_consistency_check_Gen_Mem_Inter():
    example = [
        "-GeneralParameters",
        "ProjectName SM_Ah_GKB_L",
        "SimulationTimeInPs 100000",
        "StepSizeInFs 200",
        "ReportIntervalInFs 100000",
        "Integrator LM",
        "FrictionIninvertPs 0.01",
        "setVelocitiesToTemperature true",
        "EpsilonH 4",
        "TemperatureinKelvin 300",
        "TextOutputs  XYZ PSF",
        "BinOutputs   XYZ TPK",
        "-Membrane  0",
        "MeshFile USphere_1002d_r_s0.ply",
        "NodeMass 50",
        "SpringModel N",
        "MeanCurvatureModel ItzyksonBarycentric",
        "MeanCurvatureCoeff 50",
        "Scale 1000",
        "SurfaceConstraint GWCAH 20",
        "SurfaceConstraintValue Au",
        "SurfaceConstraintCoeff 1.3",
        "VolumeConstraint N",
        "-InteractionTable",
        "M0 0",
    ]
    configurations = get_configurations(example)
    assert configurations_consistency_check(configurations) == True


def test_configurations_consistency_check_Gen_Mem():
    example = [
        "-GeneralParameters",
        "ProjectName SM_Ah_GKB_L",
        "SimulationTimeInPs 100000",
        "StepSizeInFs 200",
        "ReportIntervalInFs 100000",
        "Integrator LM",
        "FrictionIninvertPs 0.01",
        "setVelocitiesToTemperature true",
        "EpsilonH 4",
        "TemperatureinKelvin 300",
        "TextOutputs  XYZ PSF",
        "BinOutputs   XYZ TPK",
        "-Membrane  0",
        "MeshFile USphere_1002d_r_s0.ply",
        "NodeMass 50",
        "SpringModel N",
        "MeanCurvatureModel ItzyksonBarycentric",
        "MeanCurvatureCoeff 50",
        "Scale 1000",
        "SurfaceConstraint GWCAH 20",
        "SurfaceConstraintValue Au",
        "SurfaceConstraintCoeff 1.3",
        "VolumeConstraint N",
    ]
    configurations = get_configurations(example)
    assert configurations_consistency_check(configurations) == True