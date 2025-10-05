# IDE Setup for HD Library with vcpkg

This document explains how to configure various IDEs to use vcpkg with the hd library.

## Quick Start

The project now includes:
- **CMakePresets.json** - Modern CMake configuration (recommended)
- **.vscode/settings.json** - VS Code specific settings

## Visual Studio Code

### Method 1: Using CMake Presets (Recommended)

1. Install the **CMake Tools** extension
2. Open the Command Palette (`Ctrl+Shift+P`)
3. Run: `CMake: Select Configure Preset`
4. Choose: `Windows vcpkg Debug` or `Windows vcpkg Release`
5. Run: `CMake: Configure`
6. Run: `CMake: Build`

### Method 2: Using settings.json

The `.vscode/settings.json` file already configures the toolchain file automatically.

1. Delete the `build` directory
2. Reload VS Code window: `Ctrl+Shift+P` → `Developer: Reload Window`
3. CMake should auto-configure with vcpkg

## Visual Studio 2022

### Method 1: Using CMakePresets.json (Recommended)

1. Open the folder in Visual Studio
2. Visual Studio will automatically detect `CMakePresets.json`
3. In the toolbar, select the preset: `Windows vcpkg Debug`
4. Build → Build All

### Method 2: Using CMakeSettings.json

Create `CMakeSettings.json` in the project root:

```json
{
  "configurations": [
    {
      "name": "x64-Debug-vcpkg",
      "generator": "Ninja",
      "configurationType": "Debug",
      "buildRoot": "${projectDir}\\build",
      "installRoot": "${projectDir}\\out\\install\\${name}",
      "cmakeCommandArgs": "-DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake",
      "buildCommandArgs": "",
      "ctestCommandArgs": "",
      "inheritEnvironments": [ "msvc_x64_x64" ]
    }
  ]
}
```

## CLion

### Using Toolchain Settings

1. Go to: **File** → **Settings** → **Build, Execution, Deployment** → **CMake**
2. In **CMake options**, add:
   ```
   -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
   ```
3. Click **OK** and let CLion reconfigure

### Using CMakePresets.json

1. Go to: **File** → **Settings** → **Build, Execution, Deployment** → **CMake**
2. Enable: **Use CMake Presets**
3. Select preset: `windows-vcpkg-debug`

## Command Line (Any IDE)

If your IDE doesn't recognize the presets, you can manually configure:

```bash
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build .
```

## Troubleshooting

### IDE doesn't find fmt

**Problem**: IDE shows "Could not find fmt" error

**Solutions**:
1. **Delete build directory completely**
2. **Restart IDE** (important!)
3. **Verify vcpkg path**: Make sure `C:/vcpkg/scripts/buildsystems/vcpkg.cmake` exists
4. **Check CMake version**: Requires CMake 3.21+ for presets

### Wrong vcpkg instance

**Problem**: Multiple vcpkg installations causing conflicts

**Solution**: Update the toolchain file path in:
- `CMakePresets.json` (line with `CMAKE_TOOLCHAIN_FILE`)
- `.vscode/settings.json` (line with `CMAKE_TOOLCHAIN_FILE`)

Replace `C:/vcpkg/` with your actual vcpkg installation path.

### Cache issues

**Problem**: IDE keeps using old configuration

**Solution**:
1. Delete `build` directory
2. Delete `.vs` directory (Visual Studio)
3. Delete `cmake-build-*` directories (CLion)
4. Restart IDE completely
5. Reconfigure from scratch

## Verifying Configuration

After configuration, CMake output should show:
```
-- Detected vcpkg toolchain on Windows
-- ✓ Found vcpkg doctest
-- ✓ Found vcpkg fmt: 11.2.0
-- mdspan not found in vcpkg, using local copy at ../../include/mdspan/include
```

If you see these messages, vcpkg is working correctly!

## Platform-Specific Notes

### Windows
- Must use vcpkg toolchain file
- Presets: `windows-vcpkg-debug` or `windows-vcpkg-release`

### macOS
- Uses Homebrew (no toolchain file needed)
- Presets: `macos-debug` or `macos-release`
- Dependencies: `brew install doctest fmt`
