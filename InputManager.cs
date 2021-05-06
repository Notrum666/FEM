using System;
using OpenTK.Input;
using OpenTK;

namespace FEM
{
    public static class Input
    {
        private static KeyboardState prevKeyboarState;
        private static KeyboardState keyboardState;
        private static MouseState prevMouseState;
        private static MouseState mouseState;
        internal static void Init()
        {
            keyboardState = Keyboard.GetState();
            prevKeyboarState = keyboardState;
            mouseState = Mouse.GetState();
            prevMouseState = mouseState;
        }
        internal static void OnUpdateFrame()
        {
            prevKeyboarState = keyboardState;
            keyboardState = Keyboard.GetState();
            prevMouseState = mouseState;
            mouseState = Mouse.GetState();
        }
        public static bool IsKeyDown(Key key)
        {
            return keyboardState.IsKeyDown(key);
        }
        public static bool IsKeyUp(Key key)
        {
            return keyboardState.IsKeyUp(key);
        }
        public static bool IsKeyPressed(Key key)
        {
            return prevKeyboarState.IsKeyUp(key) && keyboardState.IsKeyDown(key);
        }
        public static bool IsKeyReleased(Key key)
        {
            return prevKeyboarState.IsKeyDown(key) && keyboardState.IsKeyUp(key);
        }
        public static Vector2 GetMouseDelta()
        {
            return new Vector2(mouseState.X - prevMouseState.X, mouseState.Y - prevMouseState.Y);
        }
        public static float GetScrollDelta()
        {
            return mouseState.Scroll.Y - prevMouseState.Scroll.Y;
        }
        public static bool IsMouseButtonDown(MouseButton button)
        {
            return mouseState.IsButtonDown(button);
        }
        public static bool IsMouseButtonUp(MouseButton button)
        {
            return mouseState.IsButtonUp(button);
        }
        public static bool IsMouseButtonPressed(MouseButton button)
        {
            return mouseState.IsButtonDown(button) && prevMouseState.IsButtonUp(button);
        }
        public static bool IsMouseButtonReleased(MouseButton button)
        {
            return mouseState.IsButtonDown(button) && prevMouseState.IsButtonDown(button);
        }
    }
}
