from clans3d.similarity.tool_type import ToolType
from clans3d.core.input_file_type import InputFileType
from clans3d.core.cli import parse_args
from clans3d.core.pipeline import ClansPipeline, PipelineConfig
from clans3d.utils.dependency_checks import verify_tool_dependencies
from clans3d.utils.log import setup_logging


def main():
    args = parse_args()

    setup_logging(verbose=args.verbose)

    verify_tool_dependencies(ToolType(args.tool))

    config = PipelineConfig(
        input_file=args.load,
        input_type=InputFileType(args.input_type),
        tool=ToolType(args.tool),
        foldseek_score=args.score,
    )
    pipeline = ClansPipeline(config)
    clans_file_path, cleaned_input_file_path = pipeline.run()

    print(f"CLANS file: {clans_file_path}")
    return clans_file_path, cleaned_input_file_path


if __name__ == "__main__":
    main()
